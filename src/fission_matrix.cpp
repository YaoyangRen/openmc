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

  // 初始化裂变矩阵
  fission_matrix_.resize(n_cells_, vector<double>(n_cells_, 0.0));
  current_batch_matrix_.resize(n_cells_, vector<double>(n_cells_, 0.0));
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

  // 累积到裂变矩阵: F[source_cell][fission_cell]
  {
    std::lock_guard<std::mutex> lock(data_mutex_);
    current_batch_matrix_[source_cell][fission_cell] += nu_fission;
    total_fissions_++;
  }
}

void FissionMatrix::start_new_batch(int batch_id)
{
  // 保存上一个batch的数据
  if (current_batch_id_ >= 0) {
    std::lock_guard<std::mutex> lock(data_mutex_);
    for (size_t i = 0; i < n_cells_; ++i) {
      for (size_t j = 0; j < n_cells_; ++j) {
        fission_matrix_[i][j] += current_batch_matrix_[i][j];
      }
      source_counts_[i] += current_batch_source_counts_[i];
    }
    n_realizations_++;
  }

  // 开始新batch
  current_batch_id_ = batch_id;

  // 清空当前batch数据
  for (auto& row : current_batch_matrix_) {
    std::fill(row.begin(), row.end(), 0.0);
  }
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
  std::cout << "  Estimated memory: "
            << (n_cells_ * n_cells_ * sizeof(double)) / (1024.0 * 1024.0)
            << " MB" << std::endl;

  // 归一化裂变矩阵
  std::cout << "\nNormalizing fission matrix..." << std::endl;
  vector<vector<double>> normalized_matrix(
    n_cells_, vector<double>(n_cells_, 0.0));

  for (size_t i = 0; i < n_cells_; ++i) {
    double source_total = source_counts_[i];
    if (source_total > 0.0) {
      for (size_t j = 0; j < n_cells_; ++j) {
        normalized_matrix[i][j] = fission_matrix_[i][j] / source_total;
      }
    }
  }

  // 写入HDF5文件
  hid_t file_id = file_open(filename, 'w');

  // 写入文件属性
  write_attribute(file_id, "filetype", "fission_matrix");
  write_attribute(file_id, "version", "1.0");
  write_attribute(file_id, "pitch", pitch_);
  write_attribute(file_id, "n_realizations", n_realizations_);
  write_attribute(
    file_id, "total_fissions", static_cast<int64_t>(total_fissions_.load()));
  write_attribute(
    file_id, "total_sources", static_cast<int64_t>(total_sources_.load()));

  // 写入网格信息
  write_dataset(file_id, "origin", origin_);
  write_dataset(file_id, "shape", shape_);

  // 将2D矩阵转换为1D以便写入HDF5
  vector<double> matrix_1d(n_cells_ * n_cells_);
  vector<double> normalized_matrix_1d(n_cells_ * n_cells_);

  for (size_t i = 0; i < n_cells_; ++i) {
    for (size_t j = 0; j < n_cells_; ++j) {
      size_t idx = i * n_cells_ + j;
      matrix_1d[idx] = fission_matrix_[i][j];
      normalized_matrix_1d[idx] = normalized_matrix[i][j];
    }
  }

  // 写入原始裂变矩阵数据集（2D）
  hsize_t dims[2] = {
    static_cast<hsize_t>(n_cells_), static_cast<hsize_t>(n_cells_)};
  hid_t dspace = H5Screate_simple(2, dims, nullptr);

  hid_t dset = H5Dcreate(file_id, "fission_matrix_raw", H5T_NATIVE_DOUBLE,
    dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(
    dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, matrix_1d.data());
  H5Dclose(dset);

  // 写入归一化裂变矩阵数据集（2D）
  dset = H5Dcreate(file_id, "fission_matrix_normalized", H5T_NATIVE_DOUBLE,
    dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
    normalized_matrix_1d.data());
  H5Dclose(dset);

  H5Sclose(dspace);

  // 写入源计数
  write_dataset(file_id, "source_counts", source_counts_);

  file_close(file_id);

  // 统计信息
  size_t non_zero_elements = 0;
  double max_element = 0.0;
  double sum_normalized = 0.0;

  for (size_t i = 0; i < n_cells_; ++i) {
    for (size_t j = 0; j < n_cells_; ++j) {
      if (normalized_matrix[i][j] > 0.0) {
        non_zero_elements++;
        max_element = std::max(max_element, normalized_matrix[i][j]);
        sum_normalized += normalized_matrix[i][j];
      }
    }
  }

  double sparsity =
    1.0 - static_cast<double>(non_zero_elements) / (n_cells_ * n_cells_);

  // 输出统计信息
  std::cout << "\nData Collection Summary:" << std::endl;
  std::cout << "  Total fissions recorded: " << total_fissions_ << std::endl;
  std::cout << "  Total sources recorded: " << total_sources_ << std::endl;
  std::cout << "  Number of batches: " << n_realizations_ << std::endl;

  std::cout << "\nMatrix Statistics:" << std::endl;
  std::cout << "  Non-zero elements: " << non_zero_elements << " / "
            << n_cells_ * n_cells_ << std::endl;
  std::cout << "  Sparsity: " << sparsity * 100.0 << "%" << std::endl;
  std::cout << "  Max normalized element: " << max_element << std::endl;
  std::cout << "  Sum of normalized matrix: " << sum_normalized << std::endl;

  std::cout << "\nOutput File: " << filename << std::endl;
  std::cout << std::string(70, '=') << std::endl;
}

} // namespace openmc
