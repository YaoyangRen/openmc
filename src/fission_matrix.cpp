#include "openmc/fission_matrix.h"

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/simulation.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <unordered_set>

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
  int max_batches, std::vector<double>, int score_start_batch, int n_materials)
  : grid_ {std::move(grid)}, max_batches_ {max_batches},
    score_start_batch_ {std::max(1, score_start_batch)}
{
  if (!grid_) {
    throw std::runtime_error("FissionMatrix requires a valid SharedMeshGrid.");
  }

  upper_bound_ = grid_->upper_bound();
  inv_pitch_ = grid_->inv_pitch();

  const size_t n_cells = grid_->n_cells();
  if (n_materials > 0) {
    n_materials_ = n_materials;
  } else {
    n_materials_ =
      std::max(1, static_cast<int>(model::materials.size()));
  }
  n_source_states_ = n_cells * static_cast<size_t>(n_materials_);
}

void FissionMatrix::record_source_birth(const Position& r,
  int64_t source_particle_id, int material_index, double source_weight)
{
  if (!is_scoring_current_batch())
    return;

  if (!std::isfinite(source_weight) || source_weight <= 0.0)
    return;

  const int64_t state = source_state_index(r, material_index);
  if (state < 0) {
    ++total_invalid_source_states_;
    return;
  }

  std::lock_guard<std::mutex> lock(data_mutex_);
  source_birth_states_[source_particle_id] = state;
  current_batch_source_counts_[state] += source_weight;
  total_sources_++;
}

void FissionMatrix::record_fission_site(const Position& r, double source_weight,
  int64_t source_particle_id, int material_index)
{
  if (!is_scoring_current_batch())
    return;

  if (!std::isfinite(source_weight) || source_weight <= 0.0)
    return;

  const int64_t child_state = source_state_index(r, material_index);
  if (child_state < 0) {
    ++total_invalid_fission_states_;
    return;
  }

  std::lock_guard<std::mutex> lock(data_mutex_);
  auto it = source_birth_states_.find(source_particle_id);
  if (it == source_birth_states_.end())
    return;

  StatePairKey key {it->second, child_state};
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
      for (const auto& [state, value] : current_batch_source_counts_) {
        source_counts_[state] += value;
      }
      n_realizations_++;
    } else {
      skipped_batches_++;
    }
  }

  current_batch_id_ = batch_id;
  current_batch_sparse_.clear();
  current_batch_source_counts_.clear();
  source_birth_states_.clear();
}

void FissionMatrix::compute_adjoint_source(
  const std::string& initial_guess, int max_iterations, double tolerance)
{
  std::cout << "\n" << std::string(70, '=') << std::endl;
  std::cout << "SPATIAL CLUTCH ADJOINT SOURCE COMPUTATION" << std::endl;
  std::cout << "  Mode: material-resolved source state I*(cell, material)"
            << std::endl;
  std::cout << "  Matrix score start batch: " << score_start_batch_
            << std::endl;
  std::cout << std::string(70, '=') << std::endl;

  if (fission_matrix_sparse_.empty()) {
    fatal_error("Fission matrix is empty. Cannot compute adjoint source.");
  }

  std::unordered_map<StatePairKey, double, StatePairKeyHash> normalized_matrix;
  normalized_matrix.reserve(fission_matrix_sparse_.size());
  std::unordered_set<int64_t> active_states;

  int normalized_source_rows = 0;
  for (const auto& [state, count] : source_counts_) {
    if (count > 0.0) {
      ++normalized_source_rows;
      active_states.insert(state);
    }
  }

  for (const auto& [key, value] : fission_matrix_sparse_) {
    const int64_t parent_state = key.parent_state;
    const int64_t child_state = key.child_state;
    auto source_it = source_counts_.find(parent_state);
    const double source_count =
      source_it != source_counts_.end() ? source_it->second : 0.0;
    if (source_count > 0.0) {
      normalized_matrix[key] = value / source_count;
      active_states.insert(parent_state);
      active_states.insert(child_state);
    }
  }

  if (normalized_matrix.empty()) {
    fatal_error("Normalized fission matrix is empty. Cannot compute adjoint "
                "source.");
  }

  if (initial_guess == "forward") {
    double total = 0.0;
    for (const auto& [state, count] : source_counts_) {
      total += count;
    }
    if (total > 0.0) {
      adjoint_source_.clear();
      for (const auto& [state, count] : source_counts_) {
        if (count > 0.0) {
          adjoint_source_[state] = count / total;
        }
      }
    } else {
      const double uniform = 1.0 / static_cast<double>(active_states.size());
      adjoint_source_.clear();
      for (int64_t state : active_states) {
        adjoint_source_[state] = uniform;
      }
    }
  } else {
    const double uniform = 1.0 / static_cast<double>(active_states.size());
    adjoint_source_.clear();
    for (int64_t state : active_states) {
      adjoint_source_[state] = uniform;
    }
  }

  double norm = 0.0;
  for (const auto& [state, value] : adjoint_source_) {
    norm += value;
  }
  if (norm > 0.0) {
    for (auto& [state, value] : adjoint_source_)
      value /= norm;
  }

  std::cout << "\nPerforming adjoint power iteration..." << std::endl;
  std::cout << "  Max iterations: " << max_iterations << std::endl;
  std::cout << "  Tolerance: " << tolerance << std::endl;
  std::cout << "  Normalized FM entries: " << normalized_matrix.size()
            << std::endl;
  std::cout << "  Source states with counts: " << normalized_source_rows
            << " / " << n_source_states_ << std::endl;
  std::cout << "  Scored inactive batches: " << n_realizations_ << std::endl;
  std::cout << "  Skipped inactive batches: " << skipped_batches_
            << std::endl;

  std::unordered_map<int64_t, double> I_new;
  adjoint_iterations_ = 0;
  adjoint_computed_ = false;
  adjoint_converged_ = false;
  adjoint_final_residual_ = 0.0;
  keff_reference_ = reference_keff();
  const double inv_keff = 1.0 / keff_reference_;

  for (int iter = 0; iter < max_iterations; ++iter) {
    I_new.clear();
    for (const auto& [key, F_cell] : normalized_matrix) {
      const int64_t parent_state = key.parent_state;
      const int64_t child_state = key.child_state;
      auto child_it = adjoint_source_.find(child_state);
      const double child_importance =
        child_it != adjoint_source_.end() ? child_it->second : 0.0;
      if (child_importance > 0.0) {
        I_new[parent_state] += F_cell * child_importance * inv_keff;
      }
    }

    double sum_new = 0.0;
    for (const auto& [state, value] : I_new) {
      sum_new += value;
    }
    if (sum_new == 0.0) {
      fatal_error("Spatial adjoint source collapsed to zero.");
    }

    double max_delta = 0.0;
    const double inv_sum = 1.0 / sum_new;
    std::unordered_set<int64_t> residual_states;
    residual_states.reserve(adjoint_source_.size() + I_new.size());
    for (const auto& [state, value] : adjoint_source_) {
      residual_states.insert(state);
    }
    for (const auto& [state, value] : I_new) {
      residual_states.insert(state);
    }
    for (int64_t state : residual_states) {
      auto new_it = I_new.find(state);
      const double new_val =
        new_it != I_new.end() ? new_it->second * inv_sum : 0.0;
      auto old_it = adjoint_source_.find(state);
      const double old_val =
        old_it != adjoint_source_.end() ? old_it->second : 0.0;
      max_delta =
        std::max(max_delta, std::abs(new_val - old_val));
    }
    adjoint_source_.clear();
    for (const auto& [state, value] : I_new) {
      const double new_val = value * inv_sum;
      if (new_val > 0.0) {
        adjoint_source_[state] = new_val;
      }
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

  double source_sum = 0.0;
  double max_adjoint = 0.0;
  for (const auto& [state, value] : adjoint_source_) {
    source_sum += value;
    max_adjoint = std::max(max_adjoint, value);
  }
  const int nonzero_count = get_adjoint_nonzero_states();

  std::cout << "\nAdjoint Source Statistics:" << std::endl;
  std::cout << "  Nonzero states: " << nonzero_count << " / "
            << n_source_states_
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
  const size_t n_source_states = static_cast<size_t>(n_source_states_);
  const size_t sparse_elements = fission_matrix_sparse_.size();
  const int source_nonzero_states = get_source_nonzero_states();
  const int child_nonzero_states = get_child_nonzero_states();
  const int source_nonzero_cells = get_source_nonzero_cells();
  const int child_nonzero_cells = get_child_nonzero_cells();
  const double n_cells_d = static_cast<double>(n_cells);
  const double n_states_d = static_cast<double>(n_source_states);
  const double source_state_coverage =
    n_source_states > 0 ?
      static_cast<double>(source_nonzero_states) / n_states_d :
      0.0;
  const double child_state_coverage =
    n_source_states > 0 ?
      static_cast<double>(child_nonzero_states) / n_states_d :
      0.0;
  const double source_cell_coverage =
    n_cells > 0 ? static_cast<double>(source_nonzero_cells) / n_cells_d : 0.0;
  const double child_cell_coverage =
    n_cells > 0 ? static_cast<double>(child_nonzero_cells) / n_cells_d : 0.0;
  const double nnz_fraction =
    n_source_states > 0 ?
      static_cast<double>(sparse_elements) / (n_states_d * n_states_d) :
      0.0;
  std::cout << "\nFission Matrix: " << sparse_elements
            << " non-zero elements (" << n_source_states << "x"
            << n_source_states
            << " parent_state x child_state grid) -> " << filename
            << std::endl;
  std::cout << "  Score start batch: " << score_start_batch_ << std::endl;
  std::cout << "  Scored inactive batches: " << n_realizations_ << std::endl;
  std::cout << "  Source state coverage: " << source_nonzero_states << " / "
            << n_source_states << std::endl;
  std::cout << "  Child state coverage: " << child_nonzero_states << " / "
            << n_source_states << std::endl;

  std::unordered_map<StatePairKey, double, StatePairKeyHash> normalized_sparse;
  normalized_sparse.reserve(fission_matrix_sparse_.size());
  for (const auto& [key, value] : fission_matrix_sparse_) {
    const int64_t parent_state = key.parent_state;
    auto source_it =
      source_counts_.find(parent_state);
    const double source_total =
      source_it != source_counts_.end() ? source_it->second : 0.0;
    if (source_total > 0.0) {
      normalized_sparse[key] = value / source_total;
    }
  }

  hid_t file_id = file_open(filename, 'w');
  write_attribute(file_id, "filetype", "fission_matrix_sparse");
  write_attribute(file_id, "version", "6.0");
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
  write_attribute(file_id, "total_invalid_source_states",
    static_cast<int64_t>(total_invalid_source_states_.load()));
  write_attribute(file_id, "total_invalid_fission_states",
    static_cast<int64_t>(total_invalid_fission_states_.load()));
  write_attribute(file_id, "n_cells", static_cast<int>(n_cells));
  write_attribute(file_id, "n_materials", n_materials_);
  write_attribute(file_id, "n_source_groups", n_materials_);
  write_attribute(
    file_id, "n_source_states", static_cast<int64_t>(n_source_states_));
  write_attribute(file_id, "nnz", static_cast<int64_t>(sparse_elements));
  write_attribute(file_id, "row_dim", static_cast<int64_t>(n_source_states_));
  write_attribute(file_id, "col_dim", static_cast<int64_t>(n_source_states_));
  write_attribute(file_id, "matrix_semantics",
    "row=(parent_cell,parent_material), col=(child_cell,child_material)");
  write_attribute(file_id, "source_state_definition", "cell_material");
  write_attribute(file_id, "state_indexing",
    "state=cell*n_materials+material_index");
  write_attribute(file_id, "source_states_with_counts",
    source_nonzero_states);
  write_attribute(file_id, "child_states_with_fission", child_nonzero_states);
  write_attribute(file_id, "source_state_coverage", source_state_coverage);
  write_attribute(file_id, "child_state_coverage", child_state_coverage);
  write_attribute(file_id, "source_cells_with_counts", source_nonzero_cells);
  write_attribute(file_id, "child_cells_with_fission", child_nonzero_cells);
  write_attribute(file_id, "source_cell_coverage", source_cell_coverage);
  write_attribute(file_id, "child_cell_coverage", child_cell_coverage);
  write_attribute(file_id, "nnz_fraction", nnz_fraction);

  write_dataset(file_id, "origin", grid_->origin());
  write_dataset(file_id, "shape", grid_->shape());
  std::array<double, 3> grid_pitch {
    grid_->pitch(), grid_->pitch(), grid_->pitch()};
  write_dataset(file_id, "grid_shape", grid_->shape());
  write_dataset(file_id, "grid_lower_left", grid_->origin());
  write_dataset(file_id, "grid_upper_right", grid_->upper_bound());
  write_dataset(file_id, "grid_pitch", grid_pitch);
  vector<int> material_ids;
  material_ids.reserve(n_materials_);
  for (int m = 0; m < n_materials_; ++m) {
    if (m < static_cast<int>(model::materials.size()) && model::materials[m]) {
      material_ids.push_back(model::materials[m]->id());
    } else {
      material_ids.push_back(m);
    }
  }
  write_dataset(file_id, "material_ids", material_ids);

  vector<int64_t> rows;
  vector<int64_t> cols;
  vector<int64_t> row_cells;
  vector<int> row_materials;
  vector<int64_t> col_cells;
  vector<int> col_materials;
  vector<double> values_raw;
  vector<double> values_normalized;
  rows.reserve(sparse_elements);
  cols.reserve(sparse_elements);
  row_cells.reserve(sparse_elements);
  row_materials.reserve(sparse_elements);
  col_cells.reserve(sparse_elements);
  col_materials.reserve(sparse_elements);
  values_raw.reserve(sparse_elements);
  values_normalized.reserve(sparse_elements);

  for (const auto& [key, value] : fission_matrix_sparse_) {
    const int64_t parent_state = key.parent_state;
    const int64_t child_state = key.child_state;
    rows.push_back(parent_state);
    cols.push_back(child_state);
    row_cells.push_back(parent_state / n_materials_);
    row_materials.push_back(static_cast<int>(parent_state % n_materials_));
    col_cells.push_back(child_state / n_materials_);
    col_materials.push_back(static_cast<int>(child_state % n_materials_));
    values_raw.push_back(value);

    auto it = normalized_sparse.find(key);
    values_normalized.push_back(
      it != normalized_sparse.end() ? it->second : 0.0);
  }

  write_dataset(file_id, "row_indices", rows);
  write_dataset(file_id, "col_indices", cols);
  write_dataset(file_id, "row_cells", row_cells);
  write_dataset(file_id, "row_material_indices", row_materials);
  write_dataset(file_id, "col_cells", col_cells);
  write_dataset(file_id, "col_material_indices", col_materials);
  write_dataset(file_id, "data_raw", values_raw);
  write_dataset(file_id, "data_normalized", values_normalized);
  vector<int64_t> source_count_indices;
  vector<double> source_count_values;
  source_count_indices.reserve(source_counts_.size());
  source_count_values.reserve(source_counts_.size());
  for (const auto& [state, count] : source_counts_) {
    if (count > 0.0) {
      source_count_indices.push_back(state);
      source_count_values.push_back(count);
    }
  }
  write_dataset(file_id, "source_count_indices", source_count_indices);
  write_dataset(file_id, "source_counts", source_count_values);

  if (adjoint_computed_) {
    double adj_sum = 0.0;
    for (const auto& [state, value] : adjoint_source_) {
      adj_sum += value;
    }
    if (adj_sum > 0.0) {
      for (auto& [state, value] : adjoint_source_) {
        value /= adj_sum;
      }
    }
    vector<int64_t> adjoint_source_indices;
    vector<double> adjoint_source_values;
    adjoint_source_indices.reserve(adjoint_source_.size());
    adjoint_source_values.reserve(adjoint_source_.size());
    for (const auto& [state, value] : adjoint_source_) {
      if (value > 0.0) {
        adjoint_source_indices.push_back(state);
        adjoint_source_values.push_back(value);
      }
    }
    write_dataset(file_id, "adjoint_source_indices", adjoint_source_indices);
    write_dataset(file_id, "adjoint_source", adjoint_source_values);
    write_attribute(file_id, "adjoint_source_storage", "sparse");
    write_attribute(file_id, "reference_keff", keff_reference_);
    write_attribute(
      file_id, "adjoint_iterations", static_cast<int>(adjoint_iterations_));
    write_attribute(file_id, "adjoint_converged",
      adjoint_converged_ ? 1 : 0);
    write_attribute(file_id, "adjoint_final_residual",
      adjoint_final_residual_);
    write_attribute(file_id, "adjoint_nonzero_cells",
      get_adjoint_nonzero_cells());
    write_attribute(file_id, "adjoint_nonzero_states",
      get_adjoint_nonzero_states());
    write_attribute(file_id, "adjoint_source_description",
      "Material-resolved spatial CLUTCH adjoint fission source "
      "I*(cell,material)");
    write_attribute(file_id, "adjoint_source_units", "normalized importance");
  }

  hid_t diagnostics = create_group(file_id, "diagnostics");
  write_dataset(diagnostics, "score_start_batch", score_start_batch_);
  write_dataset(diagnostics, "scored_inactive_batches", n_realizations_);
  write_dataset(diagnostics, "skipped_inactive_batches", skipped_batches_);
  write_dataset(diagnostics, "total_invalid_source_states",
    static_cast<int64_t>(total_invalid_source_states_.load()));
  write_dataset(diagnostics, "total_invalid_fission_states",
    static_cast<int64_t>(total_invalid_fission_states_.load()));
  write_dataset(diagnostics, "source_states_with_counts",
    source_nonzero_states);
  write_dataset(diagnostics, "child_states_with_fission",
    child_nonzero_states);
  write_dataset(diagnostics, "source_state_coverage",
    source_state_coverage);
  write_dataset(diagnostics, "child_state_coverage",
    child_state_coverage);
  write_dataset(diagnostics, "source_cells_with_counts",
    source_nonzero_cells);
  write_dataset(diagnostics, "child_cells_with_fission",
    child_nonzero_cells);
  write_dataset(diagnostics, "source_cell_coverage", source_cell_coverage);
  write_dataset(diagnostics, "child_cell_coverage", child_cell_coverage);
  write_dataset(diagnostics, "nnz_fraction", nnz_fraction);
  write_dataset(diagnostics, "adjoint_iterations", adjoint_iterations_);
  write_dataset(diagnostics, "adjoint_converged", adjoint_converged_ ? 1 : 0);
  write_dataset(diagnostics, "adjoint_final_residual",
    adjoint_final_residual_);
  write_dataset(diagnostics, "adjoint_nonzero_cells",
    get_adjoint_nonzero_cells());
  write_dataset(diagnostics, "adjoint_nonzero_states",
    get_adjoint_nonzero_states());
  H5Gclose(diagnostics);

  file_close(file_id);
}

int FissionMatrix::get_adjoint_nonzero_cells() const
{
  std::unordered_set<int64_t> adjoint_cells;
  for (const auto& [state, value] : adjoint_source_) {
    if (value > 1e-10) {
      adjoint_cells.insert(state / n_materials_);
    }
  }
  return static_cast<int>(adjoint_cells.size());
}

int FissionMatrix::get_adjoint_nonzero_states() const
{
  return std::count_if(adjoint_source_.begin(), adjoint_source_.end(),
    [](const auto& item) { return item.second > 1e-10; });
}

int FissionMatrix::get_source_nonzero_cells() const
{
  std::unordered_set<int64_t> source_cells;
  for (const auto& [state, count] : source_counts_) {
    if (count > 0.0) {
      source_cells.insert(state / n_materials_);
    }
  }
  return static_cast<int>(source_cells.size());
}

int FissionMatrix::get_child_nonzero_cells() const
{
  std::unordered_set<int64_t> child_cells;
  for (const auto& [key, value] : fission_matrix_sparse_) {
    if (value <= 0.0)
      continue;
    const int64_t child_state = key.child_state;
    child_cells.insert(child_state / n_materials_);
  }

  return static_cast<int>(child_cells.size());
}

int FissionMatrix::get_source_nonzero_states() const
{
  return std::count_if(source_counts_.begin(), source_counts_.end(),
    [](const auto& item) { return item.second > 0.0; });
}

int FissionMatrix::get_child_nonzero_states() const
{
  std::unordered_set<int64_t> child_states;
  for (const auto& [key, value] : fission_matrix_sparse_) {
    if (value <= 0.0)
      continue;
    child_states.insert(key.child_state);
  }

  return static_cast<int>(child_states.size());
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

int64_t FissionMatrix::source_state_index(
  const Position& r, int material_index) const
{
  if (material_index < 0 || material_index >= n_materials_)
    return -1;

  const int cell = position_to_index(r);
  if (cell < 0)
    return -1;

  return static_cast<int64_t>(cell) * n_materials_ + material_index;
}

std::unordered_map<int64_t, double> FissionMatrix::get_cell_adjoint_source()
  const
{
  std::unordered_map<int64_t, double> cell_source;
  cell_source.reserve(adjoint_source_.size());
  for (const auto& [state, value] : adjoint_source_) {
    if (value <= 0.0)
      continue;
    const int64_t cell = state / n_materials_;
    if (cell >= 0 && cell < static_cast<int64_t>(grid_->n_cells())) {
      cell_source[cell] += value;
    }
  }
  return cell_source;
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
