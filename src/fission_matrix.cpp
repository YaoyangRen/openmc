#include "openmc/fission_matrix.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/simulation.h"
#include "openmc/universe.h"

namespace openmc {

namespace {

double reference_keff()
{
  const double keff = simulation::keff;
  if (keff > 0.0 && std::isfinite(keff)) {
    return keff;
  }
  return 1.0;
}

} // namespace

FissionMatrix::FissionMatrix(
  std::shared_ptr<SharedMeshGrid> grid, int max_batches)
  : grid_ {std::move(grid)}, upper_bound_ {0.0, 0.0, 0.0}, inv_pitch_ {1.0},
    current_batch_id_ {-1}, max_batches_ {max_batches}, n_realizations_ {0},
    adjoint_computed_ {false}, adjoint_iterations_ {0}, keff_reference_ {1.0},
    enable_batch_adjoint_ {false}, adjoint_max_iter_per_batch_ {0},
    adjoint_tolerance_ {1.0e-6}, adjoint_start_batch_ {5}
{
  if (!grid_) {
    throw std::runtime_error("FissionMatrix requires a valid SharedMeshGrid.");
  }

  upper_bound_ = grid_->upper_bound();
  inv_pitch_ = grid_->inv_pitch();

  auto n_cells = grid_->n_cells();
  source_counts_.assign(n_cells, 0.0);
  current_batch_source_counts_.assign(n_cells, 0.0);
  adjoint_source_.assign(n_cells, 0.0);
  adjoint_source_batch_.assign(n_cells, 0.0);
  adjoint_source_accumulated_.assign(n_cells, 0.0);
  forward_source_.assign(n_cells, 0.0);

  if (n_cells > 0) {
    double uniform = 1.0 / static_cast<double>(n_cells);
    std::fill(adjoint_source_.begin(), adjoint_source_.end(), uniform);
  }
}

void FissionMatrix::record_source_birth(
  const Position& r, int64_t source_particle_id)
{
  int cell = position_to_index(r);
  if (cell < 0)
    return;

  std::lock_guard<std::mutex> lock(data_mutex_);
  source_birth_cells_[source_particle_id] = cell;
  current_batch_source_counts_[cell] += 1.0;
  total_sources_++;
}

void FissionMatrix::record_fission_event(
  const Position& r, double nu_fission, int64_t source_particle_id)
{
  if (nu_fission <= 0.0)
    return;

  int fission_cell = position_to_index(r);
  if (fission_cell < 0)
    return;

  std::lock_guard<std::mutex> lock(data_mutex_);
  auto it = source_birth_cells_.find(source_particle_id);
  if (it == source_birth_cells_.end())
    return;

  int source_cell = it->second;
  size_t key =
    static_cast<size_t>(source_cell) * grid_->n_cells() + fission_cell;
  current_batch_sparse_[key] += nu_fission;
  total_fissions_++;
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
    for (size_t i = 0; i < grid_->n_cells(); ++i) {
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
    size_t row = key / grid_->n_cells();
    double source_total = source_counts_[row];
    if (source_total > 0.0) {
      normalized_fm[key] = value / source_total;
    }
  }

  // 单次计算: I_new = (1/keff) × F^T × I*
  if (!normalized_fm.empty()) {
    vector<double> I_new(grid_->n_cells(), 0.0);

    for (const auto& [key, F_ij] : normalized_fm) {
      size_t i = key / grid_->n_cells(); // 源单元（行）
      size_t j = key % grid_->n_cells(); // 裂变单元（列）
      I_new[j] += F_ij * adjoint_source_[i];
    }

    const double keff = reference_keff();
    const double inv_keff = 1.0 / keff;
    for (auto& value : I_new) {
      value *= inv_keff;
    }

    double sum_new = std::accumulate(I_new.begin(), I_new.end(), 0.0);
    if (sum_new > 0.0) {
      const double inv_sum = 1.0 / sum_new;
      double max_delta = 0.0;
      for (size_t i = 0; i < grid_->n_cells(); ++i) {
        double new_val = I_new[i] * inv_sum;
        max_delta = std::max(max_delta, std::abs(new_val - adjoint_source_[i]));
        adjoint_source_[i] = new_val;
      }

      keff_reference_ = keff;
      adjoint_computed_ = true;
      adjoint_iterations_ = 1;

      std::cout << "\nAdjoint source computed using keff = " << std::fixed
                << std::setprecision(6) << keff << std::endl;
      std::cout << "  Max |ΔI*| = " << std::scientific << std::setprecision(2)
                << max_delta << std::endl;
    } else {
      std::cout << "\nWarning: adjoint source sum is zero after keff scaling"
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
    // 初始化伴随源为均匀分布
    double norm = static_cast<double>(grid_->n_cells());
    for (size_t i = 0; i < grid_->n_cells(); ++i) {
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
    for (size_t i = 0; i < grid_->n_cells(); ++i) {
      total_sources += source_counts_[i];
    }

    if (total_sources > 0.0) {
      for (size_t i = 0; i < grid_->n_cells(); ++i) {
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
    for (size_t i = 0; i < grid_->n_cells(); ++i) {
      adjoint_source_[i] /= norm;
    }
  }

  // 幂迭代法求解伴随源
  // I* = (1/k) F^T I*
  std::cout << "\nPerforming power iteration..." << std::endl;
  std::cout << "  Max iterations: " << max_iterations << std::endl;
  std::cout << "  Tolerance: " << tolerance << std::endl;

  vector<double> I_new(grid_->n_cells(), 0.0);
  double max_delta = 0.0;
  adjoint_iterations_ = 0;

  const double keff = reference_keff();
  const double inv_keff = 1.0 / keff;

  for (int iter = 0; iter < max_iterations; ++iter) {
    // 计算 F^T × I*
    // 由于 F[i][j] 存储�?key = i * n_cells + j
    // F^T[j][i] = F[i][j]
    // (F^T × I*)_j = Σ_i F^T[j][i] × I*_i = Σ_i F[i][j] × I*_i

    std::fill(I_new.begin(), I_new.end(), 0.0);

    for (const auto& [key, F_ij] : fission_matrix_sparse_) {
      size_t i = key / grid_->n_cells(); // 源单元（行）
      size_t j = key % grid_->n_cells(); // 裂变单元（列）

      // F^T[j][i] = F[i][j]
      // (F^T × I*)_j += F[i][j] × I*_i
      I_new[j] += F_ij * adjoint_source_[i];
    }

    for (auto& val : I_new) {
      val *= inv_keff;
    }

    double sum_new = std::accumulate(I_new.begin(), I_new.end(), 0.0);
    if (sum_new == 0.0) {
      std::cerr << "Error: adjoint source collapsed to zero at iteration "
                << iter << std::endl;
      break;
    }

    const double inv_sum = 1.0 / sum_new;
    max_delta = 0.0;
    for (size_t i = 0; i < grid_->n_cells(); ++i) {
      double new_val = I_new[i] * inv_sum;
      max_delta = std::max(max_delta, std::abs(new_val - adjoint_source_[i]));
      adjoint_source_[i] = new_val;
    }

    // 检查收敛性
    adjoint_iterations_ = iter + 1;

    // 50次迭代输出进度
    if ((iter + 1) % 50 == 0 || iter == 0) {
      std::cout << "  Iteration " << std::setw(4) << (iter + 1)
                << ": max |ΔI*| = " << std::scientific << std::setprecision(2)
                << max_delta << std::endl;
    }

    if (max_delta < tolerance) {
      std::cout << "\nConverged at iteration " << (iter + 1) << std::endl;
      std::cout << "  Reference keff = " << std::fixed << std::setprecision(8)
                << keff << std::endl;
      std::cout << "  Final max |ΔI*| = " << std::scientific
                << std::setprecision(2) << max_delta << std::endl;
      keff_reference_ = keff;
      adjoint_computed_ = true;
      break;
    }

    if (iter == max_iterations - 1) {
      std::cout << "\nWarning: Maximum iterations reached without convergence"
                << std::endl;
      std::cout << "  Reference keff = " << std::fixed << std::setprecision(8)
                << keff << std::endl;
      std::cout << "  Final max |ΔI*| = " << std::scientific
                << std::setprecision(2) << max_delta << std::endl;
      keff_reference_ = keff;
      adjoint_computed_ = true; // 仍标记为已计算
    }
  }

  // 计算伴随源的统计信息
  double max_adjoint = 0.0;
  double min_adjoint = 1.0e100;
  double sum_adjoint = 0.0;
  int nonzero_count = 0;

  for (size_t i = 0; i < grid_->n_cells(); ++i) {
    if (adjoint_source_[i] > 0.0) {
      max_adjoint = std::max(max_adjoint, adjoint_source_[i]);
      min_adjoint = std::min(min_adjoint, adjoint_source_[i]);
      sum_adjoint += adjoint_source_[i];
      nonzero_count++;
    }
  }

  std::cout << "\nAdjoint Source Statistics:" << std::endl;
  std::cout << "  Nonzero cells: " << nonzero_count << " / " << grid_->n_cells()
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

  // 归一化裂变矩阵（临时存储�?
  std::unordered_map<size_t, double> normalized_matrix;
  for (const auto& [key, value] : matrix_to_use) {
    size_t row = key / grid_->n_cells();
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

  vector<double> I_new(grid_->n_cells(), 0.0);
  double max_delta = 0.0;
  const double keff = reference_keff();
  const double inv_keff = 1.0 / keff;

  for (int iter = 0; iter < iterations; ++iter) {
    // 计算 F^T × I* (使用归一化矩阵)
    std::fill(I_new.begin(), I_new.end(), 0.0);

    for (const auto& [key, F_ij] : normalized_matrix) {
      size_t i = key / grid_->n_cells(); // 源单元（行）
      size_t j = key % grid_->n_cells(); // 裂变单元（列）

      // (F^T × I*)_j += F[i][j] × I*_i
      I_new[j] += F_ij * adjoint_source_[i];
    }

    for (auto& val : I_new) {
      val *= inv_keff;
    }

    double sum_new = std::accumulate(I_new.begin(), I_new.end(), 0.0);

    // 调试输出
    if (verbose && sum_new == 0.0) {
      // 检查I* 中非零元数
      int nonzero_I = 0;
      double sum_I = 0.0;
      for (size_t i = 0; i < grid_->n_cells(); ++i) {
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
        size_t i = key / grid_->n_cells();
        if (adjoint_source_[i] > 1e-15) {
          used_cells++;
          if (used_cells <= 3) {
            std::cerr << "    F[" << i << "][" << (key % grid_->n_cells())
                      << "]=" << F_ij << " * I*[" << i
                      << "]=" << adjoint_source_[i] << " = "
                      << (F_ij * adjoint_source_[i]) << std::endl;
          }
        }
      }
      std::cerr << "  Total matrix entries with nonzero I*: " << used_cells
                << std::endl;
    }

    if (sum_new == 0.0) {
      if (verbose) {
        std::cerr << "Error: adjoint source collapsed to zero at iteration "
                  << iter << std::endl;
      }
      break;
    }

    const double inv_sum = 1.0 / sum_new;
    max_delta = 0.0;
    for (size_t i = 0; i < grid_->n_cells(); ++i) {
      double new_val = I_new[i] * inv_sum;
      max_delta = std::max(max_delta, std::abs(new_val - adjoint_source_[i]));
      adjoint_source_[i] = new_val;
    }

    // 检查收敛性
    adjoint_iterations_++;

    // 详细输出（每50次迭代或收敛时）
    if (verbose && ((iter + 1) % 50 == 0 || iter == 0)) {
      std::cout << "  Iteration " << std::setw(4) << (iter + 1)
                << ": max |ΔI*| = " << std::scientific << std::setprecision(2)
                << max_delta << std::endl;
    }

    if (max_delta < adjoint_tolerance_) {
      if (verbose) {
        std::cout << "\nConverged at iteration " << (iter + 1) << std::endl;
        std::cout << "  Reference keff = " << std::fixed << std::setprecision(8)
                  << keff << std::endl;
        std::cout << "  Final max |ΔI*| = " << std::scientific
                  << std::setprecision(2) << max_delta << std::endl;
      }
      keff_reference_ = keff;
      adjoint_computed_ = true;
      break;
    }
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

  // 简洁输�?
  size_t sparse_elements = fission_matrix_sparse_.size();
  std::cout << "\nFission Matrix: " << sparse_elements << " non-zero elements ("
            << grid_->n_cells() << "x" << grid_->n_cells() << " grid) -> "
            << filename << std::endl;
  std::unordered_map<size_t, double> normalized_sparse;

  for (const auto& [key, value] : fission_matrix_sparse_) {
    size_t row = key / grid_->n_cells();
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
  write_attribute(file_id, "pitch", grid_->pitch());
  write_attribute(file_id, "n_realizations", n_realizations_);
  write_attribute(
    file_id, "total_fissions", static_cast<int64_t>(total_fissions_.load()));
  write_attribute(
    file_id, "total_sources", static_cast<int64_t>(total_sources_.load()));
  write_attribute(file_id, "n_cells", static_cast<int>(grid_->n_cells()));
  write_attribute(file_id, "nnz", static_cast<int64_t>(sparse_elements));

  // 写入网格信息
  write_dataset(file_id, "origin", grid_->origin());
  write_dataset(file_id, "shape", grid_->shape());

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
    size_t row = key / grid_->n_cells();
    size_t col = key % grid_->n_cells();
    rows.push_back(static_cast<int>(row));
    cols.push_back(static_cast<int>(col));
    values_raw.push_back(value);

    // 添加归一化
    auto it = normalized_sparse.find(key);
    values_normalized.push_back(
      it != normalized_sparse.end() ? it->second : 0.0);
  }

  // 写入稀疏矩阵数数据(COO格式)
  write_dataset(file_id, "row_indices", rows);
  write_dataset(file_id, "col_indices", cols);
  write_dataset(file_id, "data_raw", values_raw);
  write_dataset(file_id, "data_normalized", values_normalized);

  // 写入源计�?
  write_dataset(file_id, "source_counts", source_counts_);

  // 写入伴随源分布（如果已计算）
  if (adjoint_computed_) {
    // 归一化伴随源
    double adj_sum =
      std::accumulate(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
    if (adj_sum > 0.0) {
      for (auto& val : adjoint_source_) {
        val /= adj_sum;
      }
      std::cout << "Adjoint source normalized (sum = 1.0)" << std::endl;
      std::cout << "\n" << std::string(70, '=') << std::endl;
    }

    // 写入伴随源分布
    write_dataset(file_id, "adjoint_source", adjoint_source_);
    write_attribute(file_id, "reference_keff", keff_reference_);
    write_attribute(
      file_id, "adjoint_iterations", static_cast<int>(adjoint_iterations_));
    write_attribute(file_id, "adjoint_converged", adjoint_computed_);

    // 添加说明信息
    write_attribute(file_id, "adjoint_source_description",
      "Adjoint source computed from accumulated FM at end of inactive batches");
    write_attribute(file_id, "adjoint_source_units", "normalized importance");
  }

  file_close(file_id);
}

int FissionMatrix::position_to_index(const Position& r) const
{
  const auto& origin = grid_->origin();
  const auto& shape = grid_->shape();

  int indices[3];
  for (int axis = 0; axis < 3; ++axis) {
    double coord = r[axis];
    if (coord < origin[axis] || coord >= upper_bound_[axis]) {
      return -1;
    }
    int idx = static_cast<int>((coord - origin[axis]) * inv_pitch_);
    if (idx < 0) {
      idx = 0;
    } else if (idx >= shape[axis]) {
      idx = shape[axis] - 1;
    }
    indices[axis] = idx;
  }

  return (indices[0] * shape[1] + indices[1]) * shape[2] + indices[2];
}

} // namespace openmc
