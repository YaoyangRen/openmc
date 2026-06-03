#include "openmc/fission_matrix.h"

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/simulation.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

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
  int max_batches, std::vector<double>, int score_start_batch)
  : grid_ {std::move(grid)}, max_batches_ {max_batches},
    score_start_batch_ {std::max(1, score_start_batch)}
{
  if (!grid_) {
    throw std::runtime_error("FissionMatrix requires a valid SharedMeshGrid.");
  }

  upper_bound_ = grid_->upper_bound();
  inv_pitch_ = grid_->inv_pitch();

  const size_t n_cells = grid_->n_cells();
  source_counts_.assign(n_cells, 0.0);
  current_batch_source_counts_.assign(n_cells, 0.0);
  adjoint_source_.assign(n_cells, 0.0);

  if (n_cells > 0) {
    const double uniform = 1.0 / static_cast<double>(n_cells);
    std::fill(adjoint_source_.begin(), adjoint_source_.end(), uniform);
  }
}

void FissionMatrix::record_source_birth(const Position& r,
  int64_t source_particle_id, double, int, double source_weight)
{
  if (!is_scoring_current_batch())
    return;

  if (!std::isfinite(source_weight) || source_weight <= 0.0)
    return;

  const int cell = position_to_index(r);
  if (cell < 0)
    return;

  std::lock_guard<std::mutex> lock(data_mutex_);
  source_birth_cells_[source_particle_id] = cell;
  current_batch_source_counts_[cell] += source_weight;
  total_sources_++;
}

void FissionMatrix::record_fission_site(const Position& r, double source_weight,
  int64_t source_particle_id, double, int)
{
  if (!is_scoring_current_batch())
    return;

  if (!std::isfinite(source_weight) || source_weight <= 0.0)
    return;

  const int child_cell = position_to_index(r);
  if (child_cell < 0)
    return;

  std::lock_guard<std::mutex> lock(data_mutex_);
  auto it = source_birth_cells_.find(source_particle_id);
  if (it == source_birth_cells_.end())
    return;

  const size_t n_cells = grid_->n_cells();
  const size_t parent_cell = static_cast<size_t>(it->second);
  const size_t key = parent_cell * n_cells + static_cast<size_t>(child_cell);
  current_batch_sparse_[key] += source_weight;
  total_fissions_++;
}

void FissionMatrix::start_new_batch(int batch_id)
{
  if (current_batch_id_ >= 0) {
    std::lock_guard<std::mutex> lock(data_mutex_);

    if (should_score_batch(current_batch_id_)) {
      for (const auto& [key, value] : current_batch_sparse_) {
        fission_matrix_sparse_[key] += value;
      }
      for (size_t cell = 0; cell < source_counts_.size(); ++cell) {
        source_counts_[cell] += current_batch_source_counts_[cell];
      }
      n_realizations_++;
    } else {
      skipped_batches_++;
    }
  }

  current_batch_id_ = batch_id;
  current_batch_sparse_.clear();
  std::fill(current_batch_source_counts_.begin(),
    current_batch_source_counts_.end(), 0.0);
  source_birth_cells_.clear();
}

void FissionMatrix::compute_adjoint_source(
  const std::string& initial_guess, int max_iterations, double tolerance)
{
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "SPATIAL CLUTCH ADJOINT SOURCE COMPUTATION" << std::endl;
  std::cout << "  Mode: spatial source state I*(cell)" << std::endl;
  std::cout << "  Matrix score start batch: " << score_start_batch_
            << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  if (fission_matrix_sparse_.empty()) {
    fatal_error("Fission matrix is empty. Cannot compute adjoint source.");
  }

  const size_t n_cells = grid_->n_cells();
  std::unordered_map<size_t, double> normalized_matrix;
  normalized_matrix.reserve(fission_matrix_sparse_.size());

  int normalized_source_rows = 0;
  for (double count : source_counts_) {
    if (count > 0.0)
      ++normalized_source_rows;
  }

  for (const auto& [key, value] : fission_matrix_sparse_) {
    const size_t parent_cell = key / n_cells;
    const double source_count = source_counts_[parent_cell];
    if (source_count > 0.0) {
      normalized_matrix[key] = value / source_count;
    }
  }

  if (normalized_matrix.empty()) {
    fatal_error("Normalized fission matrix is empty. Cannot compute adjoint "
                "source.");
  }

  if (initial_guess == "forward") {
    const double total =
      std::accumulate(source_counts_.begin(), source_counts_.end(), 0.0);
    if (total > 0.0) {
      for (size_t cell = 0; cell < n_cells; ++cell) {
        adjoint_source_[cell] = source_counts_[cell] / total;
      }
    } else {
      const double uniform = 1.0 / static_cast<double>(n_cells);
      std::fill(adjoint_source_.begin(), adjoint_source_.end(), uniform);
    }
  } else {
    const double uniform = 1.0 / static_cast<double>(n_cells);
    std::fill(adjoint_source_.begin(), adjoint_source_.end(), uniform);
  }

  double norm =
    std::accumulate(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
  if (norm > 0.0) {
    for (auto& value : adjoint_source_)
      value /= norm;
  }

  std::cout << "\nPerforming adjoint power iteration..." << std::endl;
  std::cout << "  Max iterations: " << max_iterations << std::endl;
  std::cout << "  Tolerance: " << tolerance << std::endl;
  std::cout << "  Normalized FM entries: " << normalized_matrix.size()
            << std::endl;
  std::cout << "  Source cells with counts: " << normalized_source_rows
            << " / " << n_cells << std::endl;
  std::cout << "  Scored inactive batches: " << n_realizations_ << std::endl;
  std::cout << "  Skipped inactive batches: " << skipped_batches_
            << std::endl;

  vector<double> I_new(n_cells, 0.0);
  adjoint_iterations_ = 0;
  adjoint_computed_ = false;
  adjoint_converged_ = false;
  adjoint_final_residual_ = 0.0;
  keff_reference_ = reference_keff();
  const double inv_keff = 1.0 / keff_reference_;

  for (int iter = 0; iter < max_iterations; ++iter) {
    std::fill(I_new.begin(), I_new.end(), 0.0);
    for (const auto& [key, F_cell] : normalized_matrix) {
      const size_t parent_cell = key / n_cells;
      const size_t child_cell = key % n_cells;
      I_new[parent_cell] += F_cell * adjoint_source_[child_cell] * inv_keff;
    }

    const double sum_new = std::accumulate(I_new.begin(), I_new.end(), 0.0);
    if (sum_new == 0.0) {
      fatal_error("Spatial adjoint source collapsed to zero.");
    }

    double max_delta = 0.0;
    const double inv_sum = 1.0 / sum_new;
    for (size_t cell = 0; cell < n_cells; ++cell) {
      const double new_val = I_new[cell] * inv_sum;
      max_delta =
        std::max(max_delta, std::abs(new_val - adjoint_source_[cell]));
      adjoint_source_[cell] = new_val;
    }
    adjoint_iterations_ = iter + 1;
    adjoint_final_residual_ = max_delta;

    if ((iter + 1) % 50 == 0 || iter == 0) {
      std::cout << "  Iteration " << std::setw(4) << (iter + 1)
                << ": max |dI*| = " << std::scientific << std::setprecision(2)
                << max_delta << std::endl;
    }

    if (max_delta < tolerance) {
      std::cout << "\nConverged at iteration " << (iter + 1) << std::endl;
      adjoint_converged_ = true;
      adjoint_computed_ = true;
      break;
    }
    if (iter == max_iterations - 1) {
      std::cout << "\nWarning: Maximum iterations reached without convergence"
                << std::endl;
      adjoint_computed_ = true;
    }
  }

  const double source_sum =
    std::accumulate(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
  const double max_adjoint =
    *std::max_element(adjoint_source_.begin(), adjoint_source_.end());
  const int nonzero_count = get_adjoint_nonzero_cells();

  std::cout << "\nAdjoint Source Statistics:" << std::endl;
  std::cout << "  Nonzero cells: " << nonzero_count << " / " << n_cells
            << std::endl;
  std::cout << "  Final max |dI*|: " << std::scientific
            << std::setprecision(6) << adjoint_final_residual_ << std::endl;
  std::cout << "  Max value: " << std::scientific << std::setprecision(6)
            << max_adjoint << std::endl;
  std::cout << "  Sum: " << source_sum << std::endl;
  std::cout << std::string(70, '=') << std::endl;
}

void FissionMatrix::finalize(const std::string& filename)
{
  start_new_batch(-1);

  if (n_realizations_ == 0) {
    warning("No fission matrix data to write.");
    return;
  }

  const size_t n_cells = grid_->n_cells();
  const size_t sparse_elements = fission_matrix_sparse_.size();
  const int source_nonzero = get_source_nonzero_cells();
  const int child_nonzero = get_child_nonzero_cells();
  const double n_cells_d = static_cast<double>(n_cells);
  const double source_coverage =
    n_cells > 0 ? static_cast<double>(source_nonzero) / n_cells_d : 0.0;
  const double child_coverage =
    n_cells > 0 ? static_cast<double>(child_nonzero) / n_cells_d : 0.0;
  const double nnz_fraction =
    n_cells > 0 ? static_cast<double>(sparse_elements) / (n_cells_d * n_cells_d)
                : 0.0;
  std::cout << "\nFission Matrix: " << sparse_elements
            << " non-zero elements (" << n_cells << "x" << n_cells
            << " parent_cell x child_cell grid) -> " << filename << std::endl;
  std::cout << "  Score start batch: " << score_start_batch_ << std::endl;
  std::cout << "  Scored inactive batches: " << n_realizations_ << std::endl;
  std::cout << "  Source cell coverage: " << source_nonzero << " / " << n_cells
            << std::endl;
  std::cout << "  Child cell coverage: " << child_nonzero << " / " << n_cells
            << std::endl;

  std::unordered_map<size_t, double> normalized_sparse;
  normalized_sparse.reserve(fission_matrix_sparse_.size());
  for (const auto& [key, value] : fission_matrix_sparse_) {
    const size_t parent_cell = key / n_cells;
    const double source_total = source_counts_[parent_cell];
    if (source_total > 0.0) {
      normalized_sparse[key] = value / source_total;
    }
  }

  hid_t file_id = file_open(filename, 'w');
  write_attribute(file_id, "filetype", "fission_matrix_sparse");
  write_attribute(file_id, "version", "5.0");
  write_attribute(file_id, "storage_format", "COO");
  write_attribute(file_id, "pitch", grid_->pitch());
  write_attribute(file_id, "n_realizations", n_realizations_);
  write_attribute(file_id, "score_start_batch", score_start_batch_);
  write_attribute(file_id, "scored_inactive_batches", n_realizations_);
  write_attribute(file_id, "skipped_inactive_batches", skipped_batches_);
  write_attribute(
    file_id, "total_fissions", static_cast<int64_t>(total_fissions_.load()));
  write_attribute(
    file_id, "total_sources", static_cast<int64_t>(total_sources_.load()));
  write_attribute(file_id, "n_cells", static_cast<int>(n_cells));
  write_attribute(file_id, "n_source_groups", 1);
  write_attribute(file_id, "nnz", static_cast<int64_t>(sparse_elements));
  write_attribute(file_id, "row_dim", static_cast<int>(n_cells));
  write_attribute(file_id, "col_dim", static_cast<int>(n_cells));
  write_attribute(file_id, "matrix_semantics",
    "row=parent_cell, col=child_cell");
  write_attribute(file_id, "source_cells_with_counts", source_nonzero);
  write_attribute(file_id, "child_cells_with_fission", child_nonzero);
  write_attribute(file_id, "source_cell_coverage", source_coverage);
  write_attribute(file_id, "child_cell_coverage", child_coverage);
  write_attribute(file_id, "nnz_fraction", nnz_fraction);

  write_dataset(file_id, "origin", grid_->origin());
  write_dataset(file_id, "shape", grid_->shape());
  std::array<double, 3> grid_pitch {
    grid_->pitch(), grid_->pitch(), grid_->pitch()};
  write_dataset(file_id, "grid_shape", grid_->shape());
  write_dataset(file_id, "grid_lower_left", grid_->origin());
  write_dataset(file_id, "grid_upper_right", grid_->upper_bound());
  write_dataset(file_id, "grid_pitch", grid_pitch);

  vector<int> rows;
  vector<int> cols;
  vector<double> values_raw;
  vector<double> values_normalized;
  rows.reserve(sparse_elements);
  cols.reserve(sparse_elements);
  values_raw.reserve(sparse_elements);
  values_normalized.reserve(sparse_elements);

  for (const auto& [key, value] : fission_matrix_sparse_) {
    const size_t parent_cell = key / n_cells;
    const size_t child_cell = key % n_cells;
    rows.push_back(static_cast<int>(parent_cell));
    cols.push_back(static_cast<int>(child_cell));
    values_raw.push_back(value);

    auto it = normalized_sparse.find(key);
    values_normalized.push_back(
      it != normalized_sparse.end() ? it->second : 0.0);
  }

  write_dataset(file_id, "row_indices", rows);
  write_dataset(file_id, "col_indices", cols);
  write_dataset(file_id, "data_raw", values_raw);
  write_dataset(file_id, "data_normalized", values_normalized);
  write_dataset(file_id, "source_counts", source_counts_);

  if (adjoint_computed_) {
    double adj_sum =
      std::accumulate(adjoint_source_.begin(), adjoint_source_.end(), 0.0);
    if (adj_sum > 0.0) {
      for (auto& value : adjoint_source_) {
        value /= adj_sum;
      }
    }
    write_dataset(file_id, "adjoint_source", adjoint_source_);
    write_attribute(file_id, "reference_keff", keff_reference_);
    write_attribute(
      file_id, "adjoint_iterations", static_cast<int>(adjoint_iterations_));
    write_attribute(file_id, "adjoint_converged",
      adjoint_converged_ ? 1 : 0);
    write_attribute(file_id, "adjoint_final_residual",
      adjoint_final_residual_);
    write_attribute(file_id, "adjoint_nonzero_cells",
      get_adjoint_nonzero_cells());
    write_attribute(file_id, "adjoint_source_description",
      "Spatial CLUTCH adjoint fission source I*(cell)");
    write_attribute(file_id, "adjoint_source_units", "normalized importance");
  }

  hid_t diagnostics = create_group(file_id, "diagnostics");
  write_dataset(diagnostics, "score_start_batch", score_start_batch_);
  write_dataset(diagnostics, "scored_inactive_batches", n_realizations_);
  write_dataset(diagnostics, "skipped_inactive_batches", skipped_batches_);
  write_dataset(diagnostics, "source_cells_with_counts", source_nonzero);
  write_dataset(diagnostics, "child_cells_with_fission", child_nonzero);
  write_dataset(diagnostics, "source_cell_coverage", source_coverage);
  write_dataset(diagnostics, "child_cell_coverage", child_coverage);
  write_dataset(diagnostics, "nnz_fraction", nnz_fraction);
  write_dataset(diagnostics, "adjoint_iterations", adjoint_iterations_);
  write_dataset(diagnostics, "adjoint_converged", adjoint_converged_ ? 1 : 0);
  write_dataset(diagnostics, "adjoint_final_residual",
    adjoint_final_residual_);
  write_dataset(diagnostics, "adjoint_nonzero_cells",
    get_adjoint_nonzero_cells());
  file_close(diagnostics);

  file_close(file_id);
}

int FissionMatrix::get_adjoint_nonzero_cells() const
{
  return std::count_if(adjoint_source_.begin(), adjoint_source_.end(),
    [](double x) { return x > 1e-10; });
}

int FissionMatrix::get_source_nonzero_cells() const
{
  return std::count_if(source_counts_.begin(), source_counts_.end(),
    [](double x) { return x > 0.0; });
}

int FissionMatrix::get_child_nonzero_cells() const
{
  const size_t n_cells = grid_->n_cells();
  vector<char> child_has_fission(n_cells, 0);
  for (const auto& [key, value] : fission_matrix_sparse_) {
    if (value <= 0.0)
      continue;
    const size_t child_cell = key % n_cells;
    child_has_fission[child_cell] = 1;
  }

  return std::count(child_has_fission.begin(), child_has_fission.end(), 1);
}

int FissionMatrix::position_to_index(const Position& r) const
{
  const auto& origin = grid_->origin();
  const auto& shape = grid_->shape();

  int indices[3];
  for (int axis = 0; axis < 3; ++axis) {
    const double coord = r[axis];
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

bool FissionMatrix::should_score_batch(int batch_id) const
{
  return batch_id >= score_start_batch_;
}

bool FissionMatrix::is_scoring_current_batch() const
{
  return should_score_batch(current_batch_id_);
}

} // namespace openmc
