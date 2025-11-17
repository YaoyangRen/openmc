#include "openmc/greenfunction_mesh.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iomanip>
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
      // 静默使用手动边界，在finalize时输出警告
      std::copy(manual_lower.begin(), manual_lower.end(), llc);
      std::copy(manual_upper.begin(), manual_upper.end(), urc);
    } else {
      // 为边界添加一些边距，确保不会丢失边缘粒子
      double margin = pitch_ * 0.1;
      for (int i = 0; i < 3; ++i) {
        llc[i] -= margin;
        urc[i] += margin;
      }
      // 边界信息将在finalize时输出
    }
  } else {
    // 使用手动指定的边界
    std::copy(manual_lower.begin(), manual_lower.end(), llc);
    std::copy(manual_upper.begin(), manual_upper.end(), urc);
  }

  origin_ = {llc[0], llc[1], llc[2]};
  upper_bound_ = {urc[0], urc[1], urc[2]};

  // 计算网格尺寸
  double dx = urc[0] - llc[0];
  double dy = urc[1] - llc[1];
  double dz = urc[2] - llc[2];

  shape_[0] = std::max(1, static_cast<int>(std::ceil(dx / pitch_)) + 1);
  shape_[1] = std::max(1, static_cast<int>(std::ceil(dy / pitch_)) + 1);
  shape_[2] = std::max(1, static_cast<int>(std::ceil(dz / pitch_)) + 1);

  spatial_size_ = static_cast<size_t>(shape_[0]) * shape_[1] * shape_[2];

  // 初始化信息将在finalize时输出

  // 初始化源计数
  source_counts_.resize(spatial_size_, 0);

  // 稀疏存储：transfer_functions_sparse_ 和 cumulative_data_sparse_ 按需增长
}

int GreenFunctionMesh::position_to_cell_index(const Position& r) const
{
  constexpr double eps = 1.0e-10;

  int ix = static_cast<int>(std::floor((r.x - origin_[0] + eps) * inv_pitch_));
  int iy = static_cast<int>(std::floor((r.y - origin_[1] + eps) * inv_pitch_));
  int iz = static_cast<int>(std::floor((r.z - origin_[2] + eps) * inv_pitch_));

  // 边界检查
  if (ix < 0 || ix >= shape_[0] || iy < 0 || iy >= shape_[1] || iz < 0 ||
      iz >= shape_[2]) {
    return -1; // 超出边界
  }

  // 计算线性索引
  return ix * shape_[1] * shape_[2] + iy * shape_[2] + iz;
}

void GreenFunctionMesh::record_source_birth(
  const Position& r, int64_t source_particle_id)
{
  // 验证源粒子ID
  if (source_particle_id < 0) {
    return;
  }

  // 计算源单元索引
  int i_source = position_to_cell_index(r);

  if (i_source < 0) {
    // 源粒子在网格边界外，静默忽略
    return;
  }

  // 记录源粒子到源单元的映射
  particle_to_source_cell_[source_particle_id] = i_source;

  // 增加该源单元的计数
  source_counts_[i_source]++;

  // 稀疏存储：不需要预分配，map 会按需增长
  // transfer_functions_sparse_[i_source] 和
  // current_batch_transfer_data_sparse_[i_source] 会在 accumulate()
  // 中首次使用时自动创建
}

void GreenFunctionMesh::accumulate(
  const Position& r, double contribution, int64_t source_particle_id)
{
  // 累积传递函数 T(r_source -> r_response):
  // contribution = nu_t = (w/k_eff) × w_ufs × (ν̄Σf/Σt)
  // 表示源位置在响应位置 r 处产生的期望裂变中子数

  // 验证源粒子ID
  if (source_particle_id < 0) {
    return;
  }

  // 验证贡献值（稀疏存储：只存储非零值）
  if (!std::isfinite(contribution) || contribution <= 0.0) {
    return;
  }

  // 查找源单元
  auto it = particle_to_source_cell_.find(source_particle_id);
  if (it == particle_to_source_cell_.end()) {
    // 这个粒子没有记录源位置，忽略
    return;
  }

  int i_source = it->second;

  // 计算响应位置的单元索引
  int j_response = position_to_cell_index(r);

  if (j_response < 0) {
    dropped_contributions_++;
    return;
  }

  // 稀疏存储：只在有贡献时才创建条目
  {
    std::lock_guard<std::mutex> lock(data_mutex_);

    // 累积到当前batch的源单元数据（稀疏）
    current_batch_transfer_data_sparse_[i_source][j_response] += contribution;

    // 累积到总的传递函数（稀疏）
    cumulative_data_sparse_[j_response] += contribution;
  }

  total_contributions_++;
}
void GreenFunctionMesh::start_new_batch(int batch_id)
{
  // 保存上一个batch的数据（稀疏）
  if (current_batch_id_ >= 0) {
    for (const auto& [i_source, response_map] :
      current_batch_transfer_data_sparse_) {
      // 累积到该源单元的总传递函数中（稀疏）
      for (const auto& [j_response, value] : response_map) {
        transfer_functions_sparse_[i_source][j_response] += value;
      }
    }
  }

  // 开始新batch
  current_batch_id_ = batch_id;
  current_batch_transfer_data_sparse_.clear();
}

const vector<double>& GreenFunctionMesh::get_source_cell_data(
  int source_cell_index) const
{
  // 从稀疏存储构建稠密向量
  static vector<double> dense_data;

  auto it = transfer_functions_sparse_.find(source_cell_index);
  if (it != transfer_functions_sparse_.end()) {
    dense_data.clear();
    dense_data.resize(spatial_size_, 0.0);

    // 填充非零值
    for (const auto& [j_response, value] : it->second) {
      if (j_response >= 0 && j_response < static_cast<int>(spatial_size_)) {
        dense_data[j_response] = value;
      }
    }
    return dense_data;
  }

  // 返回空向量
  static vector<double> empty_data;
  return empty_data;
}

const vector<int>& GreenFunctionMesh::get_source_cell_indices() const
{
  static vector<int> indices;
  indices.clear();
  indices.reserve(transfer_functions_sparse_.size());

  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    indices.push_back(i_source);
  }

  std::sort(indices.begin(), indices.end());
  return indices;
}
void GreenFunctionMesh::finalize_greenfunction_mesh(
  const int batch_id, const std::string& filename)
{
  // 保存最后一个batch的数据
  start_new_batch(-1);

  if (transfer_functions_sparse_.empty()) {
    return;
  }

  // 计算有源的单元数和总源粒子数
  int n_source_cells = transfer_functions_sparse_.size();
  int total_source_particles = 0;
  for (const auto& count : source_counts_) {
    total_source_particles += count;
  }

  // 计算稀疏性统计
  size_t total_nonzero_entries = 0;
  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    total_nonzero_entries += response_map.size();
  }

  double sparsity = 100.0 * (1.0 - static_cast<double>(total_nonzero_entries) /
                                     (n_source_cells * spatial_size_));
  size_t sparse_memory_mb =
    (total_nonzero_entries * (sizeof(int) + sizeof(double)) +
      n_source_cells * 64) /
    (1024 * 1024); // 估算
  size_t dense_memory_mb =
    (n_source_cells * spatial_size_ * sizeof(double)) / (1024 * 1024);

  // 简洁输出传递函数信息
  std::cout << "\nTransfer Function T(r_s->r): " << n_source_cells
            << " source cells, " << total_source_particles
            << " source particles, " << total_contributions_ << " contributions"
            << std::endl;
  std::cout << "  Sparse storage: " << total_nonzero_entries
            << " non-zero entries (" << std::fixed << std::setprecision(2)
            << (100.0 - sparsity) << "% density, " << sparsity << "% sparsity)"
            << std::endl;
  std::cout << "  Memory saved: " << dense_memory_mb << " MB (dense) -> "
            << sparse_memory_mb << " MB (sparse), " << std::fixed
            << std::setprecision(1)
            << (100.0 * (1.0 - static_cast<double>(sparse_memory_mb) /
                                 dense_memory_mb))
            << "% reduction" << std::endl;
  std::cout << "  Output: " << filename << std::endl;

  // 创建HDF5文件，使用指定的文件名
  hid_t file_id = file_open(filename, 'w');

  // 写入文件头部信息
  write_attribute(file_id, "filetype", "transfer_function_mesh_sparse");
  write_attribute(file_id, "version", "3.0");
  write_attribute(file_id, "description",
    "Transfer function T(r_source->r_response): Sparse storage format");
  write_attribute(file_id, "pitch", pitch_);
  write_dataset(file_id, "origin", origin_);
  write_dataset(file_id, "shape", shape_);
  write_attribute(file_id, "n_source_cells", n_source_cells);
  write_attribute(file_id, "total_source_particles", total_source_particles);
  write_attribute(file_id, "total_nonzero_entries",
    static_cast<int64_t>(total_nonzero_entries));
  write_attribute(file_id, "sparsity_percent", sparsity);

  // 写入累积的传递函数（稀疏格式）
  vector<int> cumulative_indices;
  vector<double> cumulative_values;
  cumulative_indices.reserve(cumulative_data_sparse_.size());
  cumulative_values.reserve(cumulative_data_sparse_.size());

  for (const auto& [j_response, value] : cumulative_data_sparse_) {
    cumulative_indices.push_back(j_response);
    cumulative_values.push_back(value);
  }

  hid_t cumulative_group = H5Gcreate(file_id, "cumulative_transfer_function",
    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  write_dataset(cumulative_group, "indices", cumulative_indices);
  write_dataset(cumulative_group, "values", cumulative_values);
  H5Gclose(cumulative_group);

  // 收集有源的单元索引列表（排序）
  vector<int> source_cell_indices;
  source_cell_indices.reserve(n_source_cells);
  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    source_cell_indices.push_back(i_source);
  }
  std::sort(source_cell_indices.begin(), source_cell_indices.end());
  write_dataset(file_id, "source_cell_indices", source_cell_indices);

  // 收集每个源单元的粒子计数（仅有源的单元）
  vector<int> source_counts_nonzero;
  source_counts_nonzero.reserve(n_source_cells);
  for (int i_source : source_cell_indices) {
    source_counts_nonzero.push_back(source_counts_[i_source]);
  }
  write_dataset(file_id, "source_counts", source_counts_nonzero);

  // 为每个源单元创建一个组（稀疏格式）
  hid_t transfer_group = H5Gcreate(
    file_id, "transfer_functions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // 写入每个源单元的传递函数 T(i_source -> r) - 稀疏格式
  for (const auto& [i_source, response_map] : transfer_functions_sparse_) {
    std::string cell_name = "source_cell_" + std::to_string(i_source);
    hid_t cell_group = H5Gcreate(
      transfer_group, cell_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // 收集该源单元的非零响应
    vector<int> response_indices;
    vector<double> response_values;
    response_indices.reserve(response_map.size());
    response_values.reserve(response_map.size());

    for (const auto& [j_response, value] : response_map) {
      response_indices.push_back(j_response);
      response_values.push_back(value);
    }

    write_dataset(cell_group, "indices", response_indices);
    write_dataset(cell_group, "values", response_values);
    write_attribute(
      cell_group, "n_nonzero", static_cast<int>(response_map.size()));

    H5Gclose(cell_group);
  }

  H5Gclose(transfer_group);
  file_close(file_id);

  // 清理数据
  transfer_functions_sparse_.clear();
  current_batch_transfer_data_sparse_.clear();
  cumulative_data_sparse_.clear();
  source_counts_.clear();
  particle_to_source_cell_.clear();
}

} // namespace openmc