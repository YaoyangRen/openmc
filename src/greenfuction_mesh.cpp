#include "openmc/greenfunction_mesh.h"

#include <algorithm>
#include <atomic>
#include <cstdint> // for int64_t
#include <iostream>
#include <string>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/hdf5_interface.h"
#include "openmc/simulation.h"
#include "openmc/timer.h"
#include "openmc/vector.h"

namespace openmc {

GreenFunctionMesh::GreenFunctionMesh(double resolution, int max_batches)
  : pitch_(resolution), inv_pitch_(1.0 / resolution), current_batch_id_(-1),
    max_batches_(max_batches)
{

  // 手动设定边界
  double llc[3] = {0, 0, 0};
  double urc[3] = {10, 10, 10};
  origin_ = {llc[0], llc[1], llc[2]};

  // 计算网格尺寸
  double dx = urc[0] - llc[0];
  double dy = urc[1] - llc[1];
  double dz = urc[2] - llc[2];

  // 修正：确保完全覆盖边界
  shape_[0] = std::max(1, static_cast<int>(std::ceil(dx / pitch_)) + 1);
  shape_[1] = std::max(1, static_cast<int>(std::ceil(dy / pitch_)) + 1);
  shape_[2] = std::max(1, static_cast<int>(std::ceil(dz / pitch_)) + 1);

  size_t spatial_size = static_cast<size_t>(shape_[0]) * shape_[1] * shape_[2];

  // 初始化batch数据存储
  batch_data_.reserve(max_batches_);
  current_batch_data_.resize(spatial_size, 0.0);
}

void GreenFunctionMesh::accumulate(const Position& r, double contribution)
{
  constexpr double eps = 1.0e-10;

  int ix = static_cast<int>(std::floor((r.x - origin_[0] + eps) * inv_pitch_));
  int iy = static_cast<int>(std::floor((r.y - origin_[1] + eps) * inv_pitch_));
  int iz = static_cast<int>(std::floor((r.z - origin_[2] + eps) * inv_pitch_));

  // 边界检查
  if (ix >= 0 && ix < shape_[0] && iy >= 0 && iy < shape_[1] && iz >= 0 &&
      iz < shape_[2]) {
    size_t index = static_cast<size_t>(ix) +
                   static_cast<size_t>(shape_[0]) *
                     (static_cast<size_t>(iy) + static_cast<size_t>(shape_[1]) *
                                                  static_cast<size_t>(iz));

    // 向当前batch数据累积
#pragma omp atomic
    current_batch_data_[index] += contribution;

    total_contributions_++;
  } else {
    dropped_contributions_++;
  }
}

void GreenFunctionMesh::finalize_greenfunction_mesh(const int batch_id)
{
  // std::cout << "=== FINALIZING GreenFunctionMesh ===" << std::endl;

  // 保存最后一个batch的数据
  if (!current_batch_data_.empty()) {
    double batch_total = 0.0;
    for (const auto& val : current_batch_data_) {
      batch_total += val;
    }

    if (batch_total > 0.0) {
      batch_data_.push_back(current_batch_data_);
    }
  }

  if (batch_data_.empty()) {
    return;
  }

  // 创建HDF5文件
  hid_t file_id = file_open("green_function_data.h5", 'w');

  // 写入文件头部信息
  write_attribute(file_id, "filetype", "green_function_mesh");
  write_attribute(file_id, "version", "1.0");
  write_attribute(file_id, "pitch", pitch_);
  write_dataset(file_id, "origin", origin_);
  write_dataset(file_id, "shape", shape_);
  write_attribute(file_id, "n_batches", static_cast<int>(batch_data_.size()));

  // 写入每个batch的数据
  for (size_t i = 0; i < batch_data_.size(); ++i) {
    std::string dataset_name = "batch_" + std::to_string(i);
    write_dataset(file_id, dataset_name.c_str(), batch_data_[i]);
  }

  // 计算并写入累积数据
  if (!batch_data_.empty()) {
    vector<double> cumulative_data = batch_data_[0]; // 复制第一个batch

    for (size_t i = 1; i < batch_data_.size(); ++i) {
      for (size_t j = 0; j < cumulative_data.size(); ++j) {
        cumulative_data[j] += batch_data_[i][j];
      }
    }

    write_dataset(file_id, "cumulative_data", cumulative_data);
  }

  file_close(file_id);
  // std::cout << "Green function data written to green_function_data.h5"
  //           << std::endl;

  // 清理数据
  batch_data_.clear();
  current_batch_data_.clear();
}

void GreenFunctionMesh::start_new_batch(int batch_id)
{
  // 如果有之前的batch数据，保存它
  if (current_batch_id_ >= 0 && !current_batch_data_.empty()) {
    // 检查是否有非零数据
    double batch_total = 0.0;
    for (const auto& val : current_batch_data_) {
      batch_total += val;
    }

    if (batch_total > 0.0) {
      batch_data_.push_back(current_batch_data_);
      // std::cout << "Saved batch " << current_batch_id_
      //           << " with total value: " << batch_total << std::endl;
    }
  }

  // 开始新batch
  current_batch_id_ = batch_id;
  std::fill(current_batch_data_.begin(), current_batch_data_.end(), 0.0);
}

} // namespace openmc