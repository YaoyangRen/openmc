#include "openmc/greenfunction_mesh.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/simulation.h"
#include "openmc/timer.h"
#include "openmc/universe.h"
#include "openmc/vector.h"

namespace openmc {

GreenFunctionMesh::GreenFunctionMesh(double resolution, int max_batches,
  bool auto_bounds, const std::array<double, 3>& manual_lower,
  const std::array<double, 3>& manual_upper)
  : pitch_(resolution), inv_pitch_(1.0 / resolution), current_batch_id_(-1),
    max_batches_(max_batches)
{

  double llc[3];
  double urc[3];

  if (auto_bounds) {
    // 自动从根宇宙获取边界
    auto bbox = model::universes.at(model::root_universe)->bounding_box();

    llc[0] = bbox.xmin;
    llc[1] = bbox.ymin;
    llc[2] = bbox.zmin;

    urc[0] = bbox.xmax;
    urc[1] = bbox.ymax;
    urc[2] = bbox.zmax;

    // 检查是否有无限边界
    bool has_infinite = false;
    for (int i = 0; i < 3; ++i) {
      if (llc[i] <= -INFTY || urc[i] >= INFTY) {
        has_infinite = true;
        break;
      }
    }

    if (has_infinite) {
      std::cerr << "Warning: Model has infinite boundaries. Using manual "
                   "bounds instead."
                << std::endl;
      std::copy(manual_lower.begin(), manual_lower.end(), llc);
      std::copy(manual_upper.begin(), manual_upper.end(), urc);
    } else {
      // 为边界添加一些边距，确保不会丢失边缘粒子
      double margin = pitch_ * 0.1;
      for (int i = 0; i < 3; ++i) {
        llc[i] -= margin;
        urc[i] += margin;
      }

      std::cout << "Auto-detected geometry bounds from root universe:"
                << std::endl;
      std::cout << "  Original: [" << bbox.xmin << "," << bbox.ymin << ","
                << bbox.zmin << "] to [" << bbox.xmax << "," << bbox.ymax << ","
                << bbox.zmax << "]" << std::endl;
      std::cout << "  With margin (" << margin << "): [" << llc[0] << ","
                << llc[1] << "," << llc[2] << "] to [" << urc[0] << ","
                << urc[1] << "," << urc[2] << "]" << std::endl;
    }
  } else {
    // 使用手动指定的边界
    std::copy(manual_lower.begin(), manual_lower.end(), llc);
    std::copy(manual_upper.begin(), manual_upper.end(), urc);
  }

  origin_ = {llc[0], llc[1], llc[2]};

  // 计算网格尺寸
  double dx = urc[0] - llc[0];
  double dy = urc[1] - llc[1];
  double dz = urc[2] - llc[2];

  shape_[0] = std::max(1, static_cast<int>(std::ceil(dx / pitch_)) + 1);
  shape_[1] = std::max(1, static_cast<int>(std::ceil(dy / pitch_)) + 1);
  shape_[2] = std::max(1, static_cast<int>(std::ceil(dz / pitch_)) + 1);

  spatial_size_ = static_cast<size_t>(shape_[0]) * shape_[1] * shape_[2];

  // 调试信息
  std::cout << "GreenFunctionMesh initialized:" << std::endl;
  std::cout << "  Bounds: [" << llc[0] << "," << llc[1] << "," << llc[2]
            << "] to [" << urc[0] << "," << urc[1] << "," << urc[2] << "]"
            << std::endl;
  std::cout << "  Pitch: " << pitch_ << std::endl;
  std::cout << "  Shape: [" << shape_[0] << "," << shape_[1] << "," << shape_[2]
            << "]" << std::endl;
  std::cout << "  Spatial size: " << spatial_size_ << std::endl;

  // 初始化累积数据
  cumulative_data_.resize(spatial_size_, 0.0);
}

void GreenFunctionMesh::accumulate(
  const Position& r, double contribution, int64_t source_particle_id)
{
  constexpr double eps = 1.0e-10;

  // 验证源粒子ID
  if (source_particle_id < 0) {
    return; // 静默忽略无效ID，避免过多调试输出
  }

  // 验证贡献值
  if (!std::isfinite(contribution) || contribution < 0.0) {
    return; // 忽略无效的贡献值
  }

  int ix = static_cast<int>(std::floor((r.x - origin_[0] + eps) * inv_pitch_));
  int iy = static_cast<int>(std::floor((r.y - origin_[1] + eps) * inv_pitch_));
  int iz = static_cast<int>(std::floor((r.z - origin_[2] + eps) * inv_pitch_));

  // 边界检查
  if (ix >= 0 && ix < shape_[0] && iy >= 0 && iy < shape_[1] && iz >= 0 &&
      iz < shape_[2]) {
    // 计算线性索引 EG: index = ix + shape_[0] * (iy + shape_[1] * iz);
    size_t index = static_cast<size_t>(ix) +
                   static_cast<size_t>(shape_[0]) *
                     (static_cast<size_t>(iy) + static_cast<size_t>(shape_[1]) *
                                                  static_cast<size_t>(iz));

    // 额外的安全检查，确保索引在有效范围内
    if (index >= spatial_size_) {
      // 这不应该发生，但为了安全起见
      std::cerr << "Warning: GreenFunctionMesh index out of bounds: " << index
                << " >= " << spatial_size_ << " (ix=" << ix << ", iy=" << iy
                << ", iz=" << iz << ")" << std::endl;
      return;
    }

    // 线程安全地获取或创建粒子数据
    vector<double>* particle_data_ptr = nullptr;

    {
      std::lock_guard<std::mutex> lock(data_mutex_);

      auto it = current_batch_particle_data_.find(source_particle_id);
      if (it != current_batch_particle_data_.end()) {
        particle_data_ptr = &(it->second);

        // 确保数据大小正确
        if (particle_data_ptr->size() != spatial_size_) {
          particle_data_ptr->resize(spatial_size_, 0.0);
        }
      } else {
        // 创建新的粒子数据条目
        auto result = current_batch_particle_data_.emplace(
          source_particle_id, vector<double>(spatial_size_, 0.0));
        particle_data_ptr = &(result.first->second);
      }
    }

    // 在锁外进行原子累积操作
    if (particle_data_ptr == nullptr) {
      return;
    }

    // 为特定源粒子累积贡献
#pragma omp atomic
    (*particle_data_ptr)[index] += contribution;

    // 同时累积到总的格林函数中
#pragma omp atomic
    cumulative_data_[index] += contribution;

    // total_contributions_++;
  } else {
    // dropped_contributions_++;
  }
}

void GreenFunctionMesh::start_new_batch(int batch_id)
{
  // 保存上一个batch的数据
  if (current_batch_id_ >= 0) {
    for (const auto& [particle_id, data] : current_batch_particle_data_) {
      // 检查是否有非零数据
      double particle_total = 0.0;
      for (const auto& val : data) {
        particle_total += val;
      }

      if (particle_total > 0.0) {
        // 如果这个源粒子的格林函数还不存在，创建它
        if (particle_green_functions_.find(particle_id) ==
            particle_green_functions_.end()) {
          particle_green_functions_[particle_id].resize(spatial_size_, 0.0);
        }

        // 累积到该源粒子的总格林函数中
        for (size_t i = 0; i < spatial_size_; ++i) {
          // 安全检查
          if (i >= data.size()) {
            std::cerr
              << "Warning: Data index out of bounds in start_new_batch: " << i
              << " >= " << data.size() << std::endl;
            break;
          }
          if (i >= particle_green_functions_[particle_id].size()) {
            std::cerr
              << "Warning: Particle green function index out of bounds: " << i
              << " >= " << particle_green_functions_[particle_id].size()
              << std::endl;
            break;
          }
          particle_green_functions_[particle_id][i] += data[i];
        }
      }
    }
  }

  // 开始新batch
  current_batch_id_ = batch_id;
  current_batch_particle_data_.clear();
}

const vector<double>& GreenFunctionMesh::get_particle_data(
  int64_t source_particle_id) const
{
  auto it = particle_green_functions_.find(source_particle_id);
  if (it != particle_green_functions_.end()) {
    return it->second;
  }

  // 返回空向量或抛出异常
  static vector<double> empty_data;
  return empty_data;
}

void GreenFunctionMesh::finalize_greenfunction_mesh(
  const int batch_id, const std::string& filename)
{
  // 保存最后一个batch的数据
  start_new_batch(-1); // 这会保存当前batch的数据

  if (particle_green_functions_.empty()) {
    return;
  }

  // 创建HDF5文件，使用指定的文件名
  hid_t file_id = file_open(filename, 'w');

  // 写入文件头部信息
  write_attribute(file_id, "filetype", "green_function_mesh_per_particle");
  write_attribute(file_id, "version", "1.0");
  write_attribute(file_id, "pitch", pitch_);
  write_dataset(file_id, "origin", origin_);
  write_dataset(file_id, "shape", shape_);
  write_attribute(file_id, "n_source_particles",
    static_cast<int>(particle_green_functions_.size()));

  // 写入累积的总格林函数
  write_dataset(file_id, "cumulative_green_function", cumulative_data_);

  // 为每个源粒子创建一个组
  hid_t particles_group = H5Gcreate(
    file_id, "source_particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // 写入每个源粒子的格林函数
  for (const auto& [particle_id, data] : particle_green_functions_) {
    std::string particle_name = "particle_" + std::to_string(particle_id);
    write_dataset(particles_group, particle_name.c_str(), data);
  }

  // 写入源粒子ID列表
  vector<int64_t> particle_ids;
  particle_ids.reserve(particle_green_functions_.size());
  for (const auto& [particle_id, data] : particle_green_functions_) {
    particle_ids.push_back(particle_id);
  }
  write_dataset(file_id, "source_particle_ids", particle_ids);

  H5Gclose(particles_group);
  file_close(file_id);

  // 清理数据
  particle_green_functions_.clear();
  current_batch_particle_data_.clear();
  cumulative_data_.clear();
}

} // namespace openmc