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

FissionMatrix::FissionMatrix(std::shared_ptr<SharedMeshGrid> grid,
  int max_batches, std::vector<double> energy_edges)
  : grid_ {std::move(grid)}, upper_bound_ {0.0, 0.0, 0.0}, inv_pitch_ {1.0},
    n_source_groups_ {1}, current_batch_id_ {-1}, max_batches_ {max_batches},
    n_realizations_ {0}, adjoint_computed_ {false}, adjoint_iterations_ {0},
    keff_reference_ {1.0}, enable_batch_adjoint_ {false},
    adjoint_max_iter_per_batch_ {0}, adjoint_tolerance_ {1.0e-6},
    adjoint_start_batch_ {5}
{
  if (!grid_) {
    throw std::runtime_error("FissionMatrix requires a valid SharedMeshGrid.");
  }

  upper_bound_ = grid_->upper_bound();
  inv_pitch_ = grid_->inv_pitch();

  // 处理能群边界
  if (energy_edges.size() >= 2) {
    if (!std::is_sorted(energy_edges.begin(), energy_edges.end())) {
      throw std::runtime_error(
        "FissionMatrix: energy_edges must be sorted ascending.");
    }
    source_energy_edges_ = std::move(energy_edges);
    n_source_groups_ = static_cast<int>(source_energy_edges_.size()) - 1;
  }

  auto n_cells = grid_->n_cells();
  size_t n_source_states = static_cast<size_t>(n_cells) * n_source_groups_;

  source_counts_.assign(n_source_states, 0.0);
  current_batch_source_counts_.assign(n_source_states, 0.0);
  adjoint_source_.assign(n_cells, 0.0);
  adjoint_source_grouped_.assign(n_source_states, 0.0);
  adjoint_source_batch_.assign(n_cells, 0.0);
  adjoint_source_accumulated_.assign(n_cells, 0.0);
  forward_source_.assign(n_cells, 0.0);

  if (n_cells > 0) {
    double uniform = 1.0 / static_cast<double>(n_source_states);
    std::fill(
      adjoint_source_grouped_.begin(), adjoint_source_grouped_.end(), uniform);
    double uniform_cell = 1.0 / static_cast<double>(n_cells);
    std::fill(adjoint_source_.begin(), adjoint_source_.end(), uniform_cell);
  }
}

void FissionMatrix::record_source_birth(
  const Position& r, int64_t source_particle_id, double energy, int mg_group,
  double source_weight)
{
  if (!std::isfinite(source_weight) || source_weight <= 0.0)
    return;

  int cell = position_to_index(r);
  if (cell < 0)
    return;

  int g = determine_source_group(energy, mg_group);
  int source_state = cell_group_to_state(cell, g);

  std::lock_guard<std::mutex> lock(data_mutex_);
  source_birth_states_[source_particle_id] = source_state;
  current_batch_source_counts_[source_state] += source_weight;
  total_sources_++;
}

void FissionMatrix::record_fission_site(const Position& r, double source_weight,
  int64_t source_particle_id, double energy, int mg_group)
{
  if (!std::isfinite(source_weight) || source_weight <= 0.0)
    return;

  int child_cell = position_to_index(r);
  if (child_cell < 0)
    return;
  int child_group = determine_source_group(energy, mg_group);
  int child_source_state = cell_group_to_state(child_cell, child_group);

  std::lock_guard<std::mutex> lock(data_mutex_);
  auto it = source_birth_states_.find(source_particle_id);
  if (it == source_birth_states_.end())
    return;

  int source_state = it->second;
  size_t n_source_states =
    grid_->n_cells() * static_cast<size_t>(n_source_groups_);
  // key = parent_source_state * n_source_states + child_source_state
  size_t key = static_cast<size_t>(source_state) * n_source_states +
               static_cast<size_t>(child_source_state);
  current_batch_sparse_[key] += source_weight;
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
    // 累积源状态计数
    size_t n_source_states =
      static_cast<size_t>(grid_->n_cells()) * n_source_groups_;
    for (size_t s = 0; s < n_source_states; ++s) {
      source_counts_[s] += current_batch_source_counts_[s];
    }
    n_realizations_++;
  }

  // 开始新batch
  current_batch_id_ = batch_id;

  // 清空当前batch数据
  current_batch_sparse_.clear();
  size_t n_source_states =
    static_cast<size_t>(grid_->n_cells()) * n_source_groups_;
  std::fill(current_batch_source_counts_.begin(),
    current_batch_source_counts_.end(), 0.0);

  source_birth_states_.clear();
}

void FissionMatrix::compute_adjoint_source(
  const std::string& initial_guess, int max_iterations, double tolerance)
{
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "ADJOINT SOURCE COMPUTATION" << std::endl;
  if (n_source_groups_ > 1) {
    std::cout << "  Mode: energy-resolved, " << n_source_groups_
              << " source groups" << std::endl;
  } else {
    std::cout << "  Mode: scalar (n_source_groups=1)" << std::endl;
  }
  std::cout << std::string(70, '=') << std::endl;

  if (fission_matrix_sparse_.empty()) {
    std::cerr << "Error: Fission matrix is empty. Cannot compute adjoint "
                 "source."
              << std::endl;
    return;
  }

  const size_t n_cells = grid_->n_cells();
  const size_t n_states = n_cells * static_cast<size_t>(n_source_groups_);

  // ---------------------------------------------------------------
  // 1. Normalize the source-state fission matrix:
  //    M_norm[parent_source_state -> child_source_state]
  //      = M[...] / source_counts_[parent_source_state]
  // ---------------------------------------------------------------
  std::unordered_map<size_t, double> normalized_matrix;
  normalized_matrix.reserve(fission_matrix_sparse_.size());

  int normalized_source_rows = 0;
  for (size_t s = 0; s < n_states; ++s) {
    if (source_counts_[s] > 0.0)
      normalized_source_rows++;
  }

  for (const auto& [key, value] : fission_matrix_sparse_) {
    size_t source_state = key / n_states;
    double sc = source_counts_[source_state];
    if (sc > 0.0) {
      normalized_matrix[key] = value / sc;
    }
  }

  if (normalized_matrix.empty()) {
    std::cerr << "Error: Normalized fission matrix is empty." << std::endl;
    return;
  }

  std::cout << "\nInitializing adjoint source..." << std::endl;
  std::cout << "  Initial guess: " << initial_guess << std::endl;

  // ---------------------------------------------------------------
  // 2. Initialize I*(source_state)
  // ---------------------------------------------------------------
  if (initial_guess == "uniform") {
    double uniform = 1.0 / static_cast<double>(n_states);
    std::fill(
      adjoint_source_grouped_.begin(), adjoint_source_grouped_.end(), uniform);
  } else if (initial_guess == "forward") {
    double total = 0.0;
    for (size_t s = 0; s < n_states; ++s)
      total += source_counts_[s];
    if (total > 0.0) {
      for (size_t s = 0; s < n_states; ++s)
        adjoint_source_grouped_[s] = source_counts_[s] / total;
    } else {
      double uniform = 1.0 / static_cast<double>(n_states);
      std::fill(adjoint_source_grouped_.begin(), adjoint_source_grouped_.end(),
        uniform);
    }
  } else {
    double uniform = 1.0 / static_cast<double>(n_states);
    std::fill(
      adjoint_source_grouped_.begin(), adjoint_source_grouped_.end(), uniform);
  }

  // 归一化初始向量
  {
    double norm = std::accumulate(
      adjoint_source_grouped_.begin(), adjoint_source_grouped_.end(), 0.0);
    if (norm > 0.0) {
      for (auto& v : adjoint_source_grouped_)
        v /= norm;
    }
  }

  std::cout << "\nPerforming adjoint power iteration..." << std::endl;
  std::cout << "  Max iterations: " << max_iterations << std::endl;
  std::cout << "  Tolerance: " << tolerance << std::endl;
  std::cout << "  Normalized FM entries: " << normalized_matrix.size()
            << std::endl;
  std::cout << "  Source states with counts: " << normalized_source_rows
            << " / " << n_states << std::endl;

  // ---------------------------------------------------------------
  // 3. Adjoint power iteration using the transpose operator:
  //    I*_new(parent) =
  //      (1/k) * sum_child M_norm[parent,child] * I*(child)
  // ---------------------------------------------------------------
  vector<double> I_new(n_states, 0.0);
  double max_delta = 0.0;
  adjoint_iterations_ = 0;

  const double keff = reference_keff();
  const double inv_keff = 1.0 / keff;

  for (int iter = 0; iter < max_iterations; ++iter) {
    std::fill(I_new.begin(), I_new.end(), 0.0);
    for (const auto& [key, F_sj] : normalized_matrix) {
      size_t source_state = key / n_states;
      size_t child_state = key % n_states;
      I_new[source_state] +=
        F_sj * adjoint_source_grouped_[child_state] * inv_keff;
    }

    // Step c: 归一化
    double sum_new = std::accumulate(I_new.begin(), I_new.end(), 0.0);
    if (sum_new == 0.0) {
      std::cerr << "Error: adjoint source collapsed to zero at iteration "
                << iter << std::endl;
      break;
    }
    const double inv_sum = 1.0 / sum_new;

    // Step d: 计算收敛指标并更新
    max_delta = 0.0;
    for (size_t s = 0; s < n_states; ++s) {
      double new_val = I_new[s] * inv_sum;
      max_delta =
        std::max(max_delta, std::abs(new_val - adjoint_source_grouped_[s]));
      adjoint_source_grouped_[s] = new_val;
    }
    adjoint_iterations_ = iter + 1;

    if ((iter + 1) % 50 == 0 || iter == 0) {
      std::cout << "  Iteration " << std::setw(4) << (iter + 1)
                << ": max |ΔI*| = " << std::scientific << std::setprecision(2)
                << max_delta << std::endl;
    }

    if (max_delta < tolerance) {
      std::cout << "\nConverged at iteration " << (iter + 1) << std::endl;
      keff_reference_ = keff;
      adjoint_computed_ = true;
      break;
    }
    if (iter == max_iterations - 1) {
      std::cout << "\nWarning: Maximum iterations reached without convergence"
                << std::endl;
      keff_reference_ = keff;
      adjoint_computed_ = true;
    }
  }

  // ---------------------------------------------------------------
  // 5. 从 adjoint_source_grouped_ 折叠出标量 adjoint_source_（向后兼容）
  //    adjoint_source_[cell] = Σ_g adjoint_source_grouped_[cell*n_groups+g]
  // ---------------------------------------------------------------
  std::fill(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
  for (size_t j = 0; j < n_cells; ++j) {
    double sum_g = 0.0;
    for (int g = 0; g < n_source_groups_; ++g) {
      sum_g += adjoint_source_grouped_[j * n_source_groups_ + g];
    }
    adjoint_source_[j] = sum_g;
  }

  // 统计信息
  double max_adjoint =
    *std::max_element(adjoint_source_.begin(), adjoint_source_.end());
  double sum_adjoint =
    std::accumulate(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
  int nonzero_count = std::count_if(adjoint_source_.begin(),
    adjoint_source_.end(), [](double x) { return x > 0.0; });

  std::cout << "\nAdjoint Source Statistics:" << std::endl;
  std::cout << "  Nonzero cells: " << nonzero_count << " / " << n_cells
            << std::endl;
  std::cout << "  Max value: " << std::scientific << std::setprecision(6)
            << max_adjoint << std::endl;
  std::cout << "  Sum: " << sum_adjoint << std::endl;
  if (n_source_groups_ > 1) {
    std::cout << "  Grouped I*: " << n_states << " states (" << n_source_groups_
              << " groups × " << n_cells << " cells)" << std::endl;
  }
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

  const size_t n_cells = grid_->n_cells();
  const size_t n_states = n_cells * static_cast<size_t>(n_source_groups_);

  // Normalize the source-state fission matrix.
  std::unordered_map<size_t, double> normalized_matrix;
  for (const auto& [key, value] : matrix_to_use) {
    size_t source_state = key / n_states;
    double sc = source_to_use[source_state];
    if (sc > 0.0) {
      normalized_matrix[key] = value / sc;
    }
  }

  if (normalized_matrix.empty()) {
    if (verbose) {
      std::cerr
        << "Warning: Normalized matrix is empty (all source counts are zero)"
        << std::endl;
    }
    return;
  }

  if (verbose) {
    std::cout << "         Normalized FM: " << normalized_matrix.size()
              << " entries, starting adjoint iteration..." << std::endl;
  }

  vector<double> I_new(n_states, 0.0);
  double max_delta = 0.0;
  const double keff = reference_keff();
  const double inv_keff = 1.0 / keff;

  for (int iter = 0; iter < iterations; ++iter) {
    std::fill(I_new.begin(), I_new.end(), 0.0);
    for (const auto& [key, F_sj] : normalized_matrix) {
      size_t source_state = key / n_states;
      size_t child_state = key % n_states;
      I_new[source_state] +=
        F_sj * adjoint_source_grouped_[child_state] * inv_keff;
    }

    double sum_new = std::accumulate(I_new.begin(), I_new.end(), 0.0);
    if (sum_new == 0.0) {
      if (verbose) {
        std::cerr << "Error: adjoint source collapsed to zero at iteration "
                  << iter << std::endl;
      }
      break;
    }

    const double inv_sum = 1.0 / sum_new;
    max_delta = 0.0;
    for (size_t s = 0; s < n_states; ++s) {
      double new_val = I_new[s] * inv_sum;
      max_delta =
        std::max(max_delta, std::abs(new_val - adjoint_source_grouped_[s]));
      adjoint_source_grouped_[s] = new_val;
    }
    adjoint_iterations_++;

    if (verbose && ((iter + 1) % 50 == 0 || iter == 0)) {
      std::cout << "  Iteration " << std::setw(4) << (iter + 1)
                << ": max |ΔI*| = " << std::scientific << std::setprecision(2)
                << max_delta << std::endl;
    }

    if (max_delta < adjoint_tolerance_) {
      if (verbose) {
        std::cout << "\nConverged at iteration " << (iter + 1) << std::endl;
      }
      keff_reference_ = keff;
      adjoint_computed_ = true;
      // 折叠为标量
      std::fill(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
      for (size_t j = 0; j < n_cells; ++j) {
        for (int g = 0; g < n_source_groups_; ++g)
          adjoint_source_[j] +=
            adjoint_source_grouped_[j * n_source_groups_ + g];
      }
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

  // 简洁输出（行维度已扩展为 n_source_states）
  size_t n_cells = grid_->n_cells();
  size_t n_source_states = n_cells * static_cast<size_t>(n_source_groups_);
  size_t sparse_elements = fission_matrix_sparse_.size();
  std::cout << "\nFission Matrix: " << sparse_elements << " non-zero elements ("
            << n_source_states << "x" << n_source_states
            << " parent_source_state x child_source_state grid) -> "
            << filename << std::endl;
  std::unordered_map<size_t, double> normalized_sparse;

  for (const auto& [key, value] : fission_matrix_sparse_) {
    size_t source_state = key / n_source_states;
    double source_total = source_counts_[source_state];
    if (source_total > 0.0) {
      normalized_sparse[key] = value / source_total;
    }
  }

  // 写入HDF5文件
  hid_t file_id = file_open(filename, 'w');

  // 写入文件属性
  write_attribute(file_id, "filetype", "fission_matrix_sparse");
  write_attribute(file_id, "version", "4.0");
  write_attribute(file_id, "storage_format", "COO");
  write_attribute(file_id, "pitch", grid_->pitch());
  write_attribute(file_id, "n_realizations", n_realizations_);
  write_attribute(
    file_id, "total_fissions", static_cast<int64_t>(total_fissions_.load()));
  write_attribute(
    file_id, "total_sources", static_cast<int64_t>(total_sources_.load()));
  write_attribute(file_id, "n_cells", static_cast<int>(n_cells));
  write_attribute(file_id, "nnz", static_cast<int64_t>(sparse_elements));
  // Both COO axes use source_state = cell * n_source_groups + g_source.
  write_attribute(file_id, "row_dim", static_cast<int>(n_source_states));
  write_attribute(file_id, "col_dim", static_cast<int>(n_source_states));
  write_attribute(file_id, "matrix_semantics",
    "row=parent_source_state, col=child_source_state");

  // 写入网格信息
  write_dataset(file_id, "origin", grid_->origin());
  write_dataset(file_id, "shape", grid_->shape());
  std::array<double, 3> grid_pitch {
    grid_->pitch(), grid_->pitch(), grid_->pitch()};
  write_dataset(file_id, "grid_shape", grid_->shape());
  write_dataset(file_id, "grid_lower_left", grid_->origin());
  write_dataset(file_id, "grid_upper_right", grid_->upper_bound());
  write_dataset(file_id, "grid_pitch", grid_pitch);

  // COO data: row=parent source_state, col=sampled child source_state.
  vector<int> rows;
  vector<int> cols;
  vector<double> values_raw;
  vector<double> values_normalized;

  rows.reserve(sparse_elements);
  cols.reserve(sparse_elements);
  values_raw.reserve(sparse_elements);
  values_normalized.reserve(sparse_elements);

  for (const auto& [key, value] : fission_matrix_sparse_) {
    size_t source_state = key / n_source_states;
    size_t col = key % n_source_states;
    rows.push_back(static_cast<int>(source_state));
    cols.push_back(static_cast<int>(col));
    values_raw.push_back(value);

    // 添加归一化
    auto it = normalized_sparse.find(key);
    values_normalized.push_back(
      it != normalized_sparse.end() ? it->second : 0.0);
  }

  // 写入稀疏矩阵数据(COO格式)
  // 注意：rows 现在是 source_state 索引，不再是 source_cell
  write_dataset(file_id, "row_indices", rows);
  write_dataset(file_id, "col_indices", cols);
  write_dataset(file_id, "data_raw", values_raw);
  write_dataset(file_id, "data_normalized", values_normalized);

  // 写入源状态计数（size = n_source_states = n_cells * n_source_groups）
  write_dataset(file_id, "source_counts", source_counts_);
  write_attribute(
    file_id, "n_source_groups", static_cast<int>(n_source_groups_));
  if (!source_energy_edges_.empty()) {
    write_dataset(file_id, "source_energy_edges", source_energy_edges_);
  }

  // 写入伴随源分布（如果已计算）
  if (adjoint_computed_) {
    // 归一化标量伴随源
    double adj_sum =
      std::accumulate(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
    if (adj_sum > 0.0) {
      for (auto& val : adjoint_source_)
        val /= adj_sum;
      std::cout << "Adjoint source normalized (sum = 1.0)" << std::endl;
      std::cout << "\n" << std::string(70, '=') << std::endl;
    }

    // 写入标量伴随源（向后兼容）
    write_dataset(file_id, "adjoint_source", adjoint_source_);
    write_attribute(file_id, "reference_keff", keff_reference_);
    write_attribute(
      file_id, "adjoint_iterations", static_cast<int>(adjoint_iterations_));
    write_attribute(file_id, "adjoint_converged", adjoint_computed_);
    write_attribute(file_id, "adjoint_source_description",
      "Collapsed adjoint source I*(cell) = sum_g I*(cell,g) [backward compat]");
    write_attribute(file_id, "adjoint_source_units", "normalized importance");

    // 写入分群伴随源 I*(source_state)（新格式）
    if (n_source_groups_ > 1) {
      // 归一化分群伴随源
      double gs = std::accumulate(
        adjoint_source_grouped_.begin(), adjoint_source_grouped_.end(), 0.0);
      if (gs > 0.0) {
        for (auto& v : adjoint_source_grouped_)
          v /= gs;
      }
      write_dataset(file_id, "adjoint_source_grouped", adjoint_source_grouped_);
      write_attribute(file_id, "adjoint_source_grouped_description",
        "Energy-resolved adjoint source I*(source_state) where "
        "source_state = cell * n_source_groups + g_source");
    }
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

int FissionMatrix::determine_source_group(double energy_eV, int mg_group) const
{
  if (n_source_groups_ == 1)
    return 0;

  // MG 模式：直接使用 mg_group
  if (energy_eV < 0.0 && mg_group >= 0) {
    if (mg_group < n_source_groups_)
      return mg_group;
    return n_source_groups_ - 1;
  }

  // CE 模式：二分查找
  if (energy_eV >= 0.0 && !source_energy_edges_.empty()) {
    if (energy_eV <= source_energy_edges_.front())
      return 0;
    if (energy_eV >= source_energy_edges_.back())
      return n_source_groups_ - 1;
    auto it = std::upper_bound(
      source_energy_edges_.begin(), source_energy_edges_.end(), energy_eV);
    int g =
      static_cast<int>(std::distance(source_energy_edges_.begin(), it)) - 1;
    if (g < 0)
      g = 0;
    if (g >= n_source_groups_)
      g = n_source_groups_ - 1;
    return g;
  }

  return 0;
}

} // namespace openmc
