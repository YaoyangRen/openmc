#include "openmc/fission_matrix.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/universe.h"

namespace openmc {

FissionMatrix::FissionMatrix(double resolution, int max_batches,
  bool auto_bounds, const std::array<double, 3>& manual_lower,
  const std::array<double, 3>& manual_upper)
  : pitch_(resolution), inv_pitch_(1.0 / resolution), current_batch_id_(-1),
    max_batches_(max_batches), n_realizations_(0)
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

    // 检查无限边界
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
      // 添加边距
      double margin = pitch_ * 0.1;
      for (int i = 0; i < 3; ++i) {
        llc[i] -= margin;
        urc[i] += margin;
      }
      // 边界信息将在finalize时输出
    }
  } else {
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

  n_cells_ = static_cast<size_t>(shape_[0]) * shape_[1] * shape_[2];

  // 使用稀疏矩阵存储，初始化时不需要预分配大量内存
  // fission_matrix_sparse_ 和 current_batch_sparse_ 会按需增长
  source_counts_.resize(n_cells_, 0.0);
  current_batch_source_counts_.resize(n_cells_, 0.0);

  // 初始化信息将在finalize时输出
}

int FissionMatrix::position_to_index(const Position& r) const
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
  int index = ix + shape_[0] * (iy + shape_[1] * iz);
  return index;
}

void FissionMatrix::record_source_birth(
  const Position& r, int64_t source_particle_id)
{
  int cell_index = position_to_index(r);

  if (cell_index < 0) {
    return; // 超出网格范围
  }

  std::lock_guard<std::mutex> lock(data_mutex_);

  // 记录源粒子的出生位置
  source_birth_cells_[source_particle_id] = cell_index;

  // 增加该源区域的计数
  current_batch_source_counts_[cell_index] += 1.0;
  total_sources_++;
}

void FissionMatrix::record_fission_event(
  const Position& r, double nu_fission, int64_t source_particle_id)
{
  int fission_cell = position_to_index(r);

  if (fission_cell < 0) {
    return; // 超出网格范围
  }

  // 查找该粒子的源位置
  int source_cell = -1;
  {
    std::lock_guard<std::mutex> lock(data_mutex_);
    auto it = source_birth_cells_.find(source_particle_id);
    if (it != source_birth_cells_.end()) {
      source_cell = it->second;
    }
  }

  if (source_cell < 0) {
    return; // 未找到源位置
  }

  // 累积到稀疏裂变矩阵: F[source_cell][fission_cell]
  // 使用线性索引: key = row * n_cells + col
  {
    std::lock_guard<std::mutex> lock(data_mutex_);
    size_t key = static_cast<size_t>(source_cell) * n_cells_ + fission_cell;
    current_batch_sparse_[key] += nu_fission;
    total_fissions_++;
  }
}

void FissionMatrix::start_new_batch(int batch_id)
{
  // 保存上一个batch的数据
  if (current_batch_id_ >= 0) {
    std::lock_guard<std::mutex> lock(data_mutex_);
    // 合并当前batch的稀疏数据到累积矩阵
    for (const auto& [key, value] : current_batch_sparse_) {
      fission_matrix_sparse_[key] += value;
    }
    // 累积源计数
    for (size_t i = 0; i < n_cells_; ++i) {
      source_counts_[i] += current_batch_source_counts_[i];
    }
    n_realizations_++;
  }

  // 开始新batch
  current_batch_id_ = batch_id;

  // 清空当前batch数据
  current_batch_sparse_.clear();
  std::fill(current_batch_source_counts_.begin(),
    current_batch_source_counts_.end(), 0.0);

  source_birth_cells_.clear();
}

void FissionMatrix::finalize(const std::string& filename)
{
  // 保存最后一个batch
  start_new_batch(-1);

  if (n_realizations_ == 0) {
    std::cerr << "Warning: No fission matrix data to write" << std::endl;
    return;
  }

  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "FISSION MATRIX FINALIZATION" << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  // 输出网格配置信息
  std::cout << "\nGrid Configuration:" << std::endl;
  std::cout << "  Bounds: [" << origin_[0] << "," << origin_[1] << ","
            << origin_[2] << "] to [" << upper_bound_[0] << ","
            << upper_bound_[1] << "," << upper_bound_[2] << "]" << std::endl;
  std::cout << "  Pitch: " << pitch_ << " cm" << std::endl;
  std::cout << "  Shape: [" << shape_[0] << "," << shape_[1] << "," << shape_[2]
            << "]" << std::endl;
  std::cout << "  Total cells: " << n_cells_ << std::endl;
  std::cout << "  Matrix size: " << n_cells_ << " x " << n_cells_ << " = "
            << n_cells_ * n_cells_ << " elements" << std::endl;

  // 计算稀疏存储节省的内存
  size_t sparse_elements = fission_matrix_sparse_.size();
  size_t dense_memory = n_cells_ * n_cells_ * sizeof(double);
  size_t sparse_memory =
    sparse_elements * (2 * sizeof(size_t) + sizeof(double));
  double compression_ratio =
    1.0 - static_cast<double>(sparse_elements) / (n_cells_ * n_cells_);

  std::cout << "  Non-zero elements: " << sparse_elements << std::endl;
  std::cout << "  Sparsity: " << compression_ratio * 100.0 << "%" << std::endl;
  std::cout << "  Dense storage would need: "
            << dense_memory / (1024.0 * 1024.0) << " MB" << std::endl;
  std::cout << "  Sparse storage uses: " << sparse_memory / (1024.0 * 1024.0)
            << " MB" << std::endl;
  std::cout << "  Memory saved: "
            << (dense_memory - sparse_memory) / (1024.0 * 1024.0) << " MB"
            << std::endl;

  // 归一化裂变矩阵（仅对非零元素）
  std::cout << "\nNormalizing sparse fission matrix..." << std::endl;
  std::unordered_map<size_t, double> normalized_sparse;

  for (const auto& [key, value] : fission_matrix_sparse_) {
    size_t row = key / n_cells_;
    double source_total = source_counts_[row];
    if (source_total > 0.0) {
      normalized_sparse[key] = value / source_total;
    }
  }

  // 写入HDF5文件
  hid_t file_id = file_open(filename, 'w');

  // 写入文件属性
  write_attribute(file_id, "filetype", "fission_matrix_sparse");
  write_attribute(file_id, "version", "2.0");
  write_attribute(file_id, "storage_format", "COO"); // Coordinate format
  write_attribute(file_id, "pitch", pitch_);
  write_attribute(file_id, "n_realizations", n_realizations_);
  write_attribute(
    file_id, "total_fissions", static_cast<int64_t>(total_fissions_.load()));
  write_attribute(
    file_id, "total_sources", static_cast<int64_t>(total_sources_.load()));
  write_attribute(file_id, "n_cells", static_cast<int>(n_cells_));
  write_attribute(file_id, "nnz", static_cast<int64_t>(sparse_elements));

  // 写入网格信息
  write_dataset(file_id, "origin", origin_);
  write_dataset(file_id, "shape", shape_);

  // 准备稀疏矩阵的COO格式数据: (row, col, value)
  vector<int> rows;
  vector<int> cols;
  vector<double> values_raw;
  vector<double> values_normalized;

  rows.reserve(sparse_elements);
  cols.reserve(sparse_elements);
  values_raw.reserve(sparse_elements);
  values_normalized.reserve(sparse_elements);

  for (const auto& [key, value] : fission_matrix_sparse_) {
    size_t row = key / n_cells_;
    size_t col = key % n_cells_;
    rows.push_back(static_cast<int>(row));
    cols.push_back(static_cast<int>(col));
    values_raw.push_back(value);

    // 添加归一化值
    auto it = normalized_sparse.find(key);
    values_normalized.push_back(
      it != normalized_sparse.end() ? it->second : 0.0);
  }

  // 写入稀疏矩阵数据 (COO格式)
  write_dataset(file_id, "row_indices", rows);
  write_dataset(file_id, "col_indices", cols);
  write_dataset(file_id, "data_raw", values_raw);
  write_dataset(file_id, "data_normalized", values_normalized);

  // 写入源计数
  write_dataset(file_id, "source_counts", source_counts_);

  file_close(file_id);

  // 统计信息
  double max_element = 0.0;
  double sum_normalized = 0.0;

  for (double val : values_normalized) {
    if (val > 0.0) {
      max_element = std::max(max_element, val);
      sum_normalized += val;
    }
  }

  // 输出统计信息
  std::cout << "\nData Collection Summary:" << std::endl;
  std::cout << "  Total fissions recorded: " << total_fissions_ << std::endl;
  std::cout << "  Total sources recorded: " << total_sources_ << std::endl;
  std::cout << "  Number of batches: " << n_realizations_ << std::endl;

  std::cout << "\nSparse Matrix Statistics:" << std::endl;
  std::cout << "  Non-zero elements: " << sparse_elements << " / "
            << n_cells_ * n_cells_ << std::endl;
  std::cout << "  Sparsity: " << compression_ratio * 100.0 << "%" << std::endl;
  std::cout << "  Max normalized element: " << max_element << std::endl;
  std::cout << "  Sum of normalized matrix: " << sum_normalized << std::endl;

  std::cout << "\nStorage Format: COO (Coordinate)" << std::endl;
  std::cout << "Output File: " << filename << std::endl;
  std::cout << std::string(70, '=') << std::endl;
}

} // namespace openmc
