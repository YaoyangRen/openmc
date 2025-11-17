#include "openmc/fission_matrix.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>

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
    max_batches_(max_batches), n_realizations_(0), k_adjoint_(0.0),
    k_adjoint_batch_(0.0), k_adjoint_accumulated_(0.0),
    adjoint_computed_(false), adjoint_iterations_(0),
    enable_batch_adjoint_(false), adjoint_max_iter_per_batch_(10),
    adjoint_tolerance_(1.0e-6)
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

  // 初始化伴随源和正向源分布
  adjoint_source_.resize(n_cells_, 0.0);
  adjoint_source_batch_.resize(n_cells_, 0.0);
  adjoint_source_accumulated_.resize(n_cells_, 0.0);
  forward_source_.resize(n_cells_, 0.0);

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

void FissionMatrix::compute_batch_adjoint()
{
  // 只在调用时计算一次伴随源（用于最后一个非活跃batch）
  if (!enable_batch_adjoint_ || fission_matrix_sparse_.empty()) {
    return;
  }

  std::lock_guard<std::mutex> lock(data_mutex_);

  // 归一化累积的FM矩阵 F_norm[i][j] = F[i][j] / source_counts[i]
  std::unordered_map<size_t, double> normalized_fm;
  for (const auto& [key, value] : fission_matrix_sparse_) {
    size_t row = key / n_cells_;
    double source_total = source_counts_[row];
    if (source_total > 0.0) {
      normalized_fm[key] = value / source_total;
    }
  }

  // 单次计算: I_new = F^T × I*
  if (!normalized_fm.empty()) {
    vector<double> I_new(n_cells_, 0.0);

    for (const auto& [key, F_ij] : normalized_fm) {
      size_t i = key / n_cells_; // 源单元（行）
      size_t j = key % n_cells_; // 裂变单元（列）
      // F^T × I*: (I_new)_j += F[i][j] × I*_i
      I_new[j] += F_ij * adjoint_source_[i];
    }

    // 计算 k = sum(I_new)
    k_adjoint_ = std::accumulate(I_new.begin(), I_new.end(), 0.0);

    // 归一化: I* = I_new / k （如果 k > 0）
    if (k_adjoint_ > 0.0) {
      for (size_t i = 0; i < n_cells_; ++i) {
        adjoint_source_[i] = I_new[i] / k_adjoint_;
      }

      std::cout << "\nAdjoint source computed: k = " << std::fixed
                << std::setprecision(6) << k_adjoint_ << std::endl;

      adjoint_computed_ = true;
      adjoint_iterations_ = 1;
    } else {
      std::cout << "\nWarning: k_adjoint = 0, adjoint source not computed"
                << std::endl;
    }
  }
}
void FissionMatrix::enable_batch_adjoint_iteration(
  bool enable, int iterations_per_batch, double tolerance, int start_batch)
{
  enable_batch_adjoint_ = enable;
  adjoint_max_iter_per_batch_ = iterations_per_batch;
  adjoint_tolerance_ = tolerance;
  adjoint_start_batch_ = start_batch;

  if (enable) {
    std::cout << "\nAdjoint source: will compute at end of inactive batches"
              << std::endl;

    // 初始化伴随源为均匀分布
    double norm = static_cast<double>(n_cells_);
    for (size_t i = 0; i < n_cells_; ++i) {
      adjoint_source_[i] = 1.0 / norm;
    }
  }
}

void FissionMatrix::compute_adjoint_source(
  const std::string& initial_guess, int max_iterations, double tolerance)
{
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "ADJOINT SOURCE COMPUTATION" << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  if (fission_matrix_sparse_.empty()) {
    std::cerr << "Error: Fission matrix is empty. Cannot compute adjoint "
                 "source."
              << std::endl;
    return;
  }

  // 初始化伴随源分布 I*
  std::cout << "\nInitializing adjoint source..." << std::endl;
  std::cout << "  Initial guess: " << initial_guess << std::endl;

  if (initial_guess == "uniform") {
    // 均匀分布初始化: I* = 1 (归一化后)
    std::fill(adjoint_source_.begin(), adjoint_source_.end(), 1.0);
    std::cout << "  Using uniform distribution: I*(i) = 1.0" << std::endl;
  } else if (initial_guess == "forward") {
    // 使用正向源分布初始化: I* = S
    // 正向源分布即为 source_counts_ (已归一化)
    double total_sources = 0.0;
    for (size_t i = 0; i < n_cells_; ++i) {
      total_sources += source_counts_[i];
    }

    if (total_sources > 0.0) {
      for (size_t i = 0; i < n_cells_; ++i) {
        adjoint_source_[i] = source_counts_[i] / total_sources;
      }
      std::cout << "  Using forward source distribution" << std::endl;
    } else {
      std::cerr << "Warning: Forward source is zero, using uniform instead"
                << std::endl;
      std::fill(adjoint_source_.begin(), adjoint_source_.end(), 1.0);
    }
  } else {
    std::cerr << "Warning: Unknown initial guess '" << initial_guess
              << "', using uniform" << std::endl;
    std::fill(adjoint_source_.begin(), adjoint_source_.end(), 1.0);
  }

  // 归一化初始向量
  double norm = 0.0;
  for (double val : adjoint_source_) {
    norm += val;
  }
  if (norm > 0.0) {
    for (size_t i = 0; i < n_cells_; ++i) {
      adjoint_source_[i] /= norm;
    }
  }

  // 幂迭代法求解伴随源
  // I* = (1/k) F^T I*
  std::cout << "\nPerforming power iteration..." << std::endl;
  std::cout << "  Max iterations: " << max_iterations << std::endl;
  std::cout << "  Tolerance: " << tolerance << std::endl;

  vector<double> I_new(n_cells_, 0.0);
  double k_old = 1.0;
  k_adjoint_ = 1.0;
  adjoint_iterations_ = 0;

  for (int iter = 0; iter < max_iterations; ++iter) {
    // 计算 F^T × I*
    // 由于 F[i][j] 存储为 key = i * n_cells + j
    // F^T[j][i] = F[i][j]
    // (F^T × I*)_j = Σ_i F^T[j][i] × I*_i = Σ_i F[i][j] × I*_i

    std::fill(I_new.begin(), I_new.end(), 0.0);

    for (const auto& [key, F_ij] : fission_matrix_sparse_) {
      size_t i = key / n_cells_; // 源单元（行）
      size_t j = key % n_cells_; // 裂变单元（列）

      // F^T[j][i] = F[i][j]
      // (F^T × I*)_j += F[i][j] × I*_i
      I_new[j] += F_ij * adjoint_source_[i];
    }

    // 计算特征值 k = Σ I_new
    k_adjoint_ = 0.0;
    for (double val : I_new) {
      k_adjoint_ += val;
    }

    if (k_adjoint_ == 0.0) {
      std::cerr << "Error: k = 0 at iteration " << iter << std::endl;
      break;
    }

    // 归一化: I* = I_new / k
    for (size_t i = 0; i < n_cells_; ++i) {
      adjoint_source_[i] = I_new[i] / k_adjoint_;
    }

    // 检查收敛性
    double dk = std::abs(k_adjoint_ - k_old);
    adjoint_iterations_ = iter + 1;

    // 每50次迭代输出进度
    if ((iter + 1) % 50 == 0 || iter == 0) {
      std::cout << "  Iteration " << std::setw(4) << (iter + 1)
                << ": k_adj = " << std::setw(10) << std::fixed
                << std::setprecision(6) << k_adjoint_
                << ", dk = " << std::scientific << std::setprecision(2) << dk
                << std::endl;
    }

    if (dk < tolerance) {
      std::cout << "\nConverged at iteration " << (iter + 1) << std::endl;
      std::cout << "  Final k_adjoint = " << std::fixed << std::setprecision(8)
                << k_adjoint_ << std::endl;
      std::cout << "  Final dk = " << std::scientific << std::setprecision(2)
                << dk << std::endl;
      adjoint_computed_ = true;
      break;
    }

    k_old = k_adjoint_;

    if (iter == max_iterations - 1) {
      std::cout << "\nWarning: Maximum iterations reached without convergence"
                << std::endl;
      std::cout << "  Final k_adjoint = " << std::fixed << std::setprecision(8)
                << k_adjoint_ << std::endl;
      std::cout << "  Final dk = " << std::scientific << std::setprecision(2)
                << dk << std::endl;
      adjoint_computed_ = true; // 仍标记为已计算
    }
  }

  // 计算伴随源的统计信息
  double max_adjoint = 0.0;
  double min_adjoint = 1.0e100;
  double sum_adjoint = 0.0;
  int nonzero_count = 0;

  for (size_t i = 0; i < n_cells_; ++i) {
    if (adjoint_source_[i] > 0.0) {
      max_adjoint = std::max(max_adjoint, adjoint_source_[i]);
      min_adjoint = std::min(min_adjoint, adjoint_source_[i]);
      sum_adjoint += adjoint_source_[i];
      nonzero_count++;
    }
  }

  std::cout << "\nAdjoint Source Statistics:" << std::endl;
  std::cout << "  Nonzero cells: " << nonzero_count << " / " << n_cells_
            << std::endl;
  std::cout << "  Max value: " << std::scientific << std::setprecision(6)
            << max_adjoint << std::endl;
  std::cout << "  Min value: " << min_adjoint << std::endl;
  std::cout << "  Sum: " << sum_adjoint << std::endl;

  std::cout << std::string(70, '=') << std::endl;
}

void FissionMatrix::perform_adjoint_iteration(
  int iterations, bool verbose, bool use_current_batch)
{
  // 选择使用哪个矩阵和源计数
  const auto& matrix_to_use =
    use_current_batch ? current_batch_sparse_ : fission_matrix_sparse_;
  const auto& source_to_use =
    use_current_batch ? current_batch_source_counts_ : source_counts_;

  if (matrix_to_use.empty()) {
    if (verbose) {
      std::cerr << "Warning: "
                << (use_current_batch ? "Current batch" : "Accumulated")
                << " fission matrix is empty, skipping adjoint iteration"
                << std::endl;
    }
    return;
  }

  // 归一化裂变矩阵（临时存储）
  std::unordered_map<size_t, double> normalized_matrix;
  for (const auto& [key, value] : matrix_to_use) {
    size_t row = key / n_cells_;
    double source_total = source_to_use[row];
    if (source_total > 0.0) {
      normalized_matrix[key] = value / source_total;
    }
  }

  if (normalized_matrix.empty()) {
    if (verbose) {
      std::cerr
        << "Warning: Normalized matrix is empty (all source counts are zero)"
        << std::endl;

      // 调试：检查源计数
      int zero_sources = 0;
      for (size_t i = 0; i < source_to_use.size(); ++i) {
        if (source_to_use[i] == 0.0)
          zero_sources++;
      }
      std::cerr << "  Source cells with zero count: " << zero_sources << " / "
                << source_to_use.size() << std::endl;
    }
    return;
  }

  if (verbose) {
    std::cout << "         Normalized FM: " << normalized_matrix.size()
              << " entries, starting adjoint iteration..." << std::endl;
  }

  vector<double> I_new(n_cells_, 0.0);
  double k_old = k_adjoint_;

  for (int iter = 0; iter < iterations; ++iter) {
    // 计算 F^T × I* (使用归一化矩阵)
    std::fill(I_new.begin(), I_new.end(), 0.0);

    for (const auto& [key, F_ij] : normalized_matrix) {
      size_t i = key / n_cells_; // 源单元（行）
      size_t j = key % n_cells_; // 裂变单元（列）

      // (F^T × I*)_j += F[i][j] × I*_i
      I_new[j] += F_ij * adjoint_source_[i];
    }

    // 计算特征值 k = Σ I_new
    k_adjoint_ = 0.0;
    for (double val : I_new) {
      k_adjoint_ += val;
    }

    // 调试输出
    if (verbose && k_adjoint_ == 0.0) {
      // 检查 I* 中非零元素
      int nonzero_I = 0;
      double sum_I = 0.0;
      for (size_t i = 0; i < n_cells_; ++i) {
        if (adjoint_source_[i] > 0.0) {
          nonzero_I++;
          sum_I += adjoint_source_[i];
        }
      }
      std::cerr << "  Debug at iter " << iter << ": I* has " << nonzero_I
                << " nonzero cells, sum=" << sum_I << std::endl;

      // 检查哪些单元的 I* 被用到
      int used_cells = 0;
      for (const auto& [key, F_ij] : normalized_matrix) {
        size_t i = key / n_cells_;
        if (adjoint_source_[i] > 1e-15) {
          used_cells++;
          if (used_cells <= 3) {
            std::cerr << "    F[" << i << "][" << (key % n_cells_)
                      << "]=" << F_ij << " * I*[" << i
                      << "]=" << adjoint_source_[i] << " = "
                      << (F_ij * adjoint_source_[i]) << std::endl;
          }
        }
      }
      std::cerr << "  Total matrix entries with nonzero I*: " << used_cells
                << std::endl;
    }

    if (k_adjoint_ == 0.0) {
      if (verbose) {
        std::cerr << "Error: k = 0 at iteration " << iter << std::endl;
      }
      break;
    }

    // 归一化: I* = I_new / k
    for (size_t i = 0; i < n_cells_; ++i) {
      adjoint_source_[i] = I_new[i] / k_adjoint_;
    }

    // 检查收敛性
    double dk = std::abs(k_adjoint_ - k_old);
    adjoint_iterations_++;

    // 详细输出（每50次迭代或收敛时）
    if (verbose && ((iter + 1) % 50 == 0 || iter == 0)) {
      std::cout << "  Iteration " << std::setw(4) << (iter + 1)
                << ": k_adj = " << std::setw(10) << std::fixed
                << std::setprecision(6) << k_adjoint_
                << ", dk = " << std::scientific << std::setprecision(2) << dk
                << std::endl;
    }

    if (dk < adjoint_tolerance_) {
      if (verbose) {
        std::cout << "\nConverged at iteration " << (iter + 1) << std::endl;
        std::cout << "  Final k_adjoint = " << std::fixed
                  << std::setprecision(8) << k_adjoint_ << std::endl;
        std::cout << "  Final dk = " << std::scientific << std::setprecision(2)
                  << dk << std::endl;
      }
      adjoint_computed_ = true;
      break;
    }

    k_old = k_adjoint_;
  }
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

  // 写入伴随源分布（如果已计算）
  if (adjoint_computed_) {
    // 归一化伴随源，使其总和为1
    double adj_sum =
      std::accumulate(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
    if (adj_sum > 0.0) {
      for (auto& val : adjoint_source_) {
        val /= adj_sum;
      }
      std::cout << "Adjoint source normalized (sum = 1.0)" << std::endl;
    }

    // 写入伴随源分布
    write_dataset(file_id, "adjoint_source", adjoint_source_);
    write_attribute(file_id, "k_adjoint", k_adjoint_);
    write_attribute(file_id, "adjoint_converged", adjoint_computed_);

    // 添加说明信息
    write_attribute(file_id, "adjoint_source_description",
      "Adjoint source computed from accumulated FM at end of inactive batches");
    write_attribute(file_id, "adjoint_source_units", "normalized importance");
  }

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

  // 输出伴随源信息（如果已计算）
  if (adjoint_computed_) {
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "ADJOINT SOURCE SUMMARY" << std::endl;
    std::cout << std::string(70, '=') << std::endl;

    // 统计非零单元
    int nonzero_cells = std::count_if(adjoint_source_.begin(),
      adjoint_source_.end(), [](double x) { return x > 1e-10; });

    double max_importance =
      *std::max_element(adjoint_source_.begin(), adjoint_source_.end());

    std::cout << "k_adjoint = " << std::fixed << std::setprecision(6)
              << k_adjoint_ << std::endl;
    std::cout << "Nonzero cells: " << nonzero_cells << " / " << n_cells_
              << std::endl;
    std::cout << "Max importance: " << std::scientific << std::setprecision(4)
              << max_importance << std::endl;
    std::cout << "Saved to '" << filename << "': adjoint_source" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
  }

  std::cout << "\nStorage Format: COO (Coordinate)" << std::endl;
  std::cout << "Output File: " << filename << std::endl;
  std::cout << std::string(70, '=') << std::endl;
}

} // namespace openmc
