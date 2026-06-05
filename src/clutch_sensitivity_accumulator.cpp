#include "openmc/clutch_sensitivity_accumulator.h"

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/mesh_init.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/tallies/derivative.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_set>
#include <utility>

#ifdef OPENMC_MPI
#include <mpi.h>
#endif

namespace openmc {

namespace {

const char* variable_name(DerivativeVariable variable)
{
  switch (variable) {
  case DerivativeVariable::DENSITY:
    return "density";
  case DerivativeVariable::NUCLIDE_DENSITY:
    return "nuclide_density";
  case DerivativeVariable::TEMPERATURE:
    return "temperature";
  }
  return "unknown";
}

int material_index_from_id(int material_id)
{
  auto it = model::material_map.find(material_id);
  if (it == model::material_map.end()) {
    fatal_error("CLUTCH sensitivity derivative references material " +
                std::to_string(material_id) + ", which is not defined.");
  }
  return it->second;
}

} // namespace

ClutchSensitivityAccumulator::ClutchSensitivityAccumulator(
  std::shared_ptr<SharedMeshGrid> grid, vector<int> derivative_indices,
  std::string method)
  : grid_ {std::move(grid)}, derivative_indices_ {std::move(derivative_indices)},
    method_ {std::move(method)}
{
  if (!grid_) {
    fatal_error("ClutchSensitivityAccumulator requires a valid SharedMeshGrid.");
  }
  if (derivative_indices_.empty()) {
    fatal_error("CLUTCH sensitivity requires at least one TallyDerivative.");
  }

  upper_bound_ = grid_->upper_bound();
  inv_pitch_ = grid_->inv_pitch();
  initialize_parameter_metadata();
}

ClutchSensitivityAccumulator::ClutchSensitivityAccumulator(
  std::shared_ptr<SharedMeshGrid> grid, vector<int> parameter_ids,
  vector<double> parameter_values, std::string method)
  : grid_ {std::move(grid)}, derivative_indices_ {std::move(parameter_ids)},
    derivative_ids_ {derivative_indices_},
    parameter_values_ {std::move(parameter_values)}, method_ {std::move(method)}
{
  if (!grid_) {
    fatal_error("ClutchSensitivityAccumulator requires a valid SharedMeshGrid.");
  }
  if (derivative_indices_.empty() ||
      derivative_indices_.size() != parameter_values_.size()) {
    fatal_error("Synthetic CLUTCH sensitivity accumulator metadata is invalid.");
  }

  upper_bound_ = grid_->upper_bound();
  inv_pitch_ = grid_->inv_pitch();
  derivative_variables_.assign(n_parameters(), "synthetic");
  derivative_material_ids_.assign(n_parameters(), 0);
  derivative_nuclides_.assign(n_parameters(), "");
}

void ClutchSensitivityAccumulator::initialize_parameter_metadata()
{
  derivative_ids_.clear();
  derivative_variables_.clear();
  derivative_material_ids_.clear();
  derivative_nuclides_.clear();
  parameter_values_.clear();

  for (int deriv_index : derivative_indices_) {
    if (deriv_index < 0 ||
        deriv_index >= static_cast<int>(model::tally_derivs.size())) {
      fatal_error("CLUTCH sensitivity received an invalid derivative index.");
    }

    const auto& deriv = model::tally_derivs[deriv_index];
    const int material_index = material_index_from_id(deriv.diff_material);
    const auto& material = *model::materials[material_index];

    derivative_ids_.push_back(deriv.id);
    derivative_variables_.push_back(variable_name(deriv.variable));
    derivative_material_ids_.push_back(deriv.diff_material);

    double parameter_value = 0.0;
    std::string nuclide_name;
    switch (deriv.variable) {
    case DerivativeVariable::DENSITY:
      parameter_value = material.density_gpcc_;
      break;

    case DerivativeVariable::NUCLIDE_DENSITY:
      if (deriv.diff_nuclide < 0 ||
          deriv.diff_nuclide >= static_cast<int>(data::nuclides.size())) {
        fatal_error("CLUTCH sensitivity nuclide_density derivative has an "
                    "invalid nuclide index.");
      }
      nuclide_name = data::nuclides[deriv.diff_nuclide]->name_;
      for (int i = 0; i < material.nuclide_.size(); ++i) {
        if (material.nuclide_[i] == deriv.diff_nuclide) {
          parameter_value = material.atom_density_(i);
          break;
        }
      }
      if (parameter_value <= 0.0) {
        fatal_error("CLUTCH sensitivity nuclide_density derivative references "
                    "a nuclide with zero density in material " +
                    std::to_string(deriv.diff_material) + ".");
      }
      break;

    case DerivativeVariable::TEMPERATURE:
      parameter_value = material.temperature();
      break;
    }

    derivative_nuclides_.push_back(nuclide_name);
    parameter_values_.push_back(parameter_value);
  }
}

void ClutchSensitivityAccumulator::set_adjoint_source_spatial(
  const vector<double>& adjoint_source_spatial)
{
  if (adjoint_source_spatial.size() != grid_->n_cells()) {
    fatal_error("ClutchSensitivityAccumulator: spatial adjoint source size "
                "does not match the shared kinetics mesh.");
  }

  std::unordered_map<int64_t, double> sparse_source;
  for (size_t i = 0; i < adjoint_source_spatial.size(); ++i) {
    const double value = adjoint_source_spatial[i];
    if (!std::isfinite(value) || value < 0.0) {
      fatal_error("ClutchSensitivityAccumulator: spatial adjoint source "
                  "contains invalid values.");
    }
    if (value > 0.0) {
      sparse_source[static_cast<int64_t>(i)] = value;
    }
  }

  set_adjoint_source_spatial(sparse_source);
}

void ClutchSensitivityAccumulator::set_adjoint_source_spatial(
  const std::unordered_map<int64_t, double>& adjoint_source_spatial)
{
  double sum = 0.0;
  std::unordered_map<int64_t, double> sparse_source;
  sparse_source.reserve(adjoint_source_spatial.size());
  for (const auto& [cell, value] : adjoint_source_spatial) {
    if (cell < 0 || cell >= static_cast<int64_t>(grid_->n_cells())) {
      fatal_error("ClutchSensitivityAccumulator: spatial adjoint source "
                  "contains an invalid cell index.");
    }
    if (!std::isfinite(value) || value < 0.0) {
      fatal_error("ClutchSensitivityAccumulator: spatial adjoint source "
                  "contains invalid values.");
    }
    if (value > 0.0) {
      sparse_source[cell] = value;
      sum += value;
    }
  }
  if (sum <= 0.0) {
    fatal_error("ClutchSensitivityAccumulator: spatial adjoint source is zero.");
  }

  std::lock_guard<std::mutex> lock(mutex_);
  adjoint_source_spatial_ = std::move(sparse_source);
  source_ready_ = true;
}

void ClutchSensitivityAccumulator::begin_batch(int batch_id)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!source_ready_) {
    fatal_error("ClutchSensitivityAccumulator: adjoint source is not loaded.");
  }
  if (batch_active_) {
    fatal_error("ClutchSensitivityAccumulator: previous batch was not ended.");
  }

  current_batch_ = BatchScore {};
  current_batch_.batch_id = batch_id;
  current_batch_.numerator.assign(n_parameters(), 0.0);
  current_batch_.cclutch_numerator.assign(n_parameters(), 0.0);
  current_source_cells_.clear();
  current_source_counts_.clear();
  current_cclutch_transfer_total_.clear();
  current_cclutch_transfer_numerator_.clear();
  batch_active_ = true;
}

void ClutchSensitivityAccumulator::end_batch(int batch_id)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_ || current_batch_.batch_id != batch_id) {
    fatal_error("ClutchSensitivityAccumulator: invalid end_batch call.");
  }
  fold_current_cclutch_batch();
  batches_.push_back(current_batch_);
  batch_active_ = false;
  current_source_cells_.clear();
}

void ClutchSensitivityAccumulator::record_source_birth(
  const Position& r, int64_t source_particle_id)
{
  if (source_particle_id < 0) {
    return;
  }

  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    return;
  }

  const int64_t cell = position_to_index(r);
  if (cell < 0) {
    return;
  }

  auto [it, inserted] = current_source_cells_.emplace(source_particle_id, cell);
  if (inserted) {
    ++current_source_counts_[cell];
  } else if (it->second != cell) {
    auto count_it = current_source_counts_.find(it->second);
    if (count_it != current_source_counts_.end()) {
      --count_it->second;
      if (count_it->second <= 0) {
        current_source_counts_.erase(count_it);
      }
    }
    it->second = cell;
    ++current_source_counts_[cell];
  }
}

void ClutchSensitivityAccumulator::score_fission_site(
  const Particle& p, const Position& r, double site_weight)
{
  vector<double> path_terms(n_parameters(), 0.0);
  vector<double> track_terms(n_parameters(), 0.0);
  vector<double> collision_terms(n_parameters(), 0.0);
  vector<double> direct_terms(n_parameters(), 0.0);

  for (int i = 0; i < derivative_indices_.size(); ++i) {
    const int deriv_index = derivative_indices_[i];
    path_terms[i] = p.flux_derivs(deriv_index);
    track_terms[i] = p.clutch_track_derivs(deriv_index);
    collision_terms[i] = p.clutch_collision_derivs(deriv_index);
    direct_terms[i] = nu_fission_direct_derivative(p, deriv_index);
  }

  score_contribution(
    r, site_weight, path_terms, track_terms, collision_terms, direct_terms);
}

void ClutchSensitivityAccumulator::score_cclutch_fission_event(
  const Particle& p, const Position& r, double total_contribution)
{
  if (!std::isfinite(total_contribution) || total_contribution <= 0.0) {
    return;
  }

  vector<double> response_terms(n_parameters(), 0.0);
  for (int i = 0; i < derivative_indices_.size(); ++i) {
    const int deriv_index = derivative_indices_[i];
    response_terms[i] =
      p.flux_derivs(deriv_index) + nu_fission_direct_derivative(p, deriv_index);
  }

  score_cclutch_contribution(
    r, p.source_particle_id(), total_contribution, response_terms);
}

void ClutchSensitivityAccumulator::score_cclutch_contribution(const Position& r,
  int64_t source_particle_id, double total_contribution,
  const vector<double>& response_terms)
{
  if (!std::isfinite(total_contribution) || total_contribution <= 0.0) {
    return;
  }
  if (response_terms.size() != n_parameters()) {
    fatal_error("ClutchSensitivityAccumulator: C-CLUTCH response vector size "
                "does not match the number of derivatives.");
  }

  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    fatal_error("ClutchSensitivityAccumulator: C-CLUTCH event scored outside "
                "an active batch.");
  }
  if (!source_ready_) {
    fatal_error("ClutchSensitivityAccumulator: adjoint source is not loaded.");
  }

  ++current_batch_.cclutch_events;
  ++total_cclutch_events_;

  if (position_to_index(r) < 0) {
    ++current_batch_.cclutch_dropped_events;
    ++total_cclutch_dropped_events_;
    return;
  }

  auto it = current_source_cells_.find(source_particle_id);
  if (it == current_source_cells_.end()) {
    ++current_batch_.cclutch_missing_source_events;
    ++total_cclutch_missing_source_events_;
    return;
  }

  const int64_t source_cell = it->second;
  if (source_cell < 0 ||
      source_cell >= static_cast<int64_t>(grid_->n_cells())) {
    ++current_batch_.cclutch_dropped_events;
    ++total_cclutch_dropped_events_;
    return;
  }

  current_cclutch_transfer_total_[source_cell] += total_contribution;
  auto& numerator = current_cclutch_transfer_numerator_[source_cell];
  if (numerator.empty()) {
    numerator.assign(n_parameters(), 0.0);
  }
  for (size_t i = 0; i < n_parameters(); ++i) {
    numerator[i] += total_contribution * response_terms[i];
  }

  ++current_batch_.cclutch_scored_events;
  ++total_cclutch_scored_events_;
}

void ClutchSensitivityAccumulator::score_contribution(const Position& r,
  double site_weight, const vector<double>& path_terms,
  const vector<double>& track_terms, const vector<double>& collision_terms,
  const vector<double>& direct_terms)
{
  if (!std::isfinite(site_weight) || site_weight <= 0.0) {
    return;
  }

  if (path_terms.size() != n_parameters() ||
      track_terms.size() != n_parameters() ||
      collision_terms.size() != n_parameters() ||
      direct_terms.size() != n_parameters()) {
    fatal_error("ClutchSensitivityAccumulator: contribution vector size "
                "does not match the number of derivatives.");
  }

  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    fatal_error("ClutchSensitivityAccumulator: fission site scored outside an "
                "active batch.");
  }
  if (!source_ready_) {
    fatal_error("ClutchSensitivityAccumulator: adjoint source is not loaded.");
  }

  ++current_batch_.fission_sites;
  ++total_fission_sites_;

  const int cell = position_to_index(r);
  if (cell < 0) {
    ++current_batch_.dropped_sites;
    ++total_dropped_sites_;
    return;
  }

  auto importance_it = adjoint_source_spatial_.find(cell);
  if (importance_it == adjoint_source_spatial_.end() ||
      importance_it->second <= 0.0) {
    ++current_batch_.dropped_sites;
    ++total_dropped_sites_;
    return;
  }

  const double response = site_weight * importance_it->second;
  current_batch_.denominator += response;
  for (int i = 0; i < n_parameters(); ++i) {
    current_batch_.numerator[i] += response * (path_terms[i] + direct_terms[i]);
  }

  ++current_batch_.scored_sites;
  ++total_scored_sites_;
}

void ClutchSensitivityAccumulator::fold_current_cclutch_batch()
{
  const size_t n_params = n_parameters();
  std::unordered_set<int64_t> source_cells;
  source_cells.reserve(current_cclutch_transfer_total_.size() +
                       current_cclutch_transfer_numerator_.size());
  for (const auto& [source_cell, total] : current_cclutch_transfer_total_) {
    source_cells.insert(source_cell);
  }
  for (const auto& [source_cell, numerator] :
       current_cclutch_transfer_numerator_) {
    source_cells.insert(source_cell);
  }

  for (int64_t source_cell : source_cells) {
    auto total_it = current_cclutch_transfer_total_.find(source_cell);
    const double total =
      total_it != current_cclutch_transfer_total_.end() ? total_it->second :
                                                          0.0;
    bool has_numerator = false;
    auto numerator_it = current_cclutch_transfer_numerator_.find(source_cell);
    if (numerator_it != current_cclutch_transfer_numerator_.end()) {
      for (double value : numerator_it->second) {
        has_numerator = has_numerator || value != 0.0;
      }
    }
    if (total <= 0.0 && !has_numerator) {
      continue;
    }

    auto count_it = current_source_counts_.find(source_cell);
    const int source_count =
      count_it != current_source_counts_.end() ? count_it->second : 0;
    if (source_count <= 0) {
      fatal_error("ClutchSensitivityAccumulator: C-CLUTCH transfer response "
                  "exists for a source cell with zero source count.");
    }

    auto importance_it = adjoint_source_spatial_.find(source_cell);
    if (importance_it == adjoint_source_spatial_.end() ||
        importance_it->second <= 0.0) {
      continue;
    }

    const double source_weight =
      importance_it->second / static_cast<double>(source_count);
    current_batch_.cclutch_denominator += source_weight * total;
    if (numerator_it != current_cclutch_transfer_numerator_.end()) {
      if (numerator_it->second.size() != n_params) {
        fatal_error("ClutchSensitivityAccumulator: C-CLUTCH transfer numerator "
                    "has an invalid derivative vector size.");
      }
      for (size_t i = 0; i < n_params; ++i) {
        current_batch_.cclutch_numerator[i] +=
          source_weight * numerator_it->second[i];
      }
    }
  }
}

ClutchSensitivityAccumulator::MethodResult
ClutchSensitivityAccumulator::compute_result(const vector<BatchScore>& batches)
  const
{
  MethodResult result;
  const size_t n = batches.size();
  const size_t n_params = n_parameters();
  result.numerator.assign(n_params, 0.0);
  result.dlogk_dparameter.assign(n_params, 0.0);
  result.sensitivity.assign(n_params, 0.0);
  result.uncertainty.assign(n_params, 0.0);
  if (n == 0) {
    return result;
  }

  double sum_denominator = 0.0;
  for (const auto& batch : batches) {
    sum_denominator += batch.denominator;
    for (int i = 0; i < n_params; ++i) {
      result.numerator[i] += batch.numerator[i];
    }
  }
  if (sum_denominator <= 0.0) {
    return result;
  }

  const double inv_n = 1.0 / static_cast<double>(n);
  result.denominator = sum_denominator * inv_n;
  result.available = true;
  for (int i = 0; i < n_params; ++i) {
    result.numerator[i] *= inv_n;
    result.dlogk_dparameter[i] = result.numerator[i] / result.denominator;
    result.sensitivity[i] =
      parameter_values_[i] * result.dlogk_dparameter[i];
  }

  if (n > 1) {
    for (int i = 0; i < n_params; ++i) {
      const double ratio = result.dlogk_dparameter[i];
      double mean_z = 0.0;
      vector<double> z_values;
      z_values.reserve(n);
      for (const auto& batch : batches) {
        const double z = batch.numerator[i] - ratio * batch.denominator;
        z_values.push_back(z);
        mean_z += z;
      }
      mean_z *= inv_n;

      double s2 = 0.0;
      for (double z : z_values) {
        const double dz = z - mean_z;
        s2 += dz * dz;
      }
      s2 /= static_cast<double>(n - 1);
      const double sigma_ratio =
        std::sqrt(s2 / (static_cast<double>(n) * result.denominator *
                         result.denominator));
      result.uncertainty[i] = std::abs(parameter_values_[i]) * sigma_ratio;
    }
  }

  return result;
}

ClutchSensitivityAccumulator::MethodResult
ClutchSensitivityAccumulator::compute_cclutch_result(
  const vector<BatchScore>& batches) const
{
  MethodResult result;
  const size_t n = batches.size();
  const size_t n_params = n_parameters();
  result.numerator.assign(n_params, 0.0);
  result.dlogk_dparameter.assign(n_params, 0.0);
  result.sensitivity.assign(n_params, 0.0);
  result.uncertainty.assign(n_params, 0.0);
  if (n == 0) {
    return result;
  }

  double sum_denominator = 0.0;
  for (const auto& batch : batches) {
    sum_denominator += batch.cclutch_denominator;
    for (int i = 0; i < n_params; ++i) {
      result.numerator[i] += batch.cclutch_numerator[i];
    }
  }
  if (sum_denominator <= 0.0) {
    return result;
  }

  const double inv_n = 1.0 / static_cast<double>(n);
  result.denominator = sum_denominator * inv_n;
  result.available = true;
  for (int i = 0; i < n_params; ++i) {
    result.numerator[i] *= inv_n;
    result.dlogk_dparameter[i] = result.numerator[i] / result.denominator;
    result.sensitivity[i] =
      parameter_values_[i] * result.dlogk_dparameter[i];
  }

  if (n > 1) {
    for (int i = 0; i < n_params; ++i) {
      const double ratio = result.dlogk_dparameter[i];
      double mean_z = 0.0;
      vector<double> z_values;
      z_values.reserve(n);
      for (const auto& batch : batches) {
        const double z =
          batch.cclutch_numerator[i] - ratio * batch.cclutch_denominator;
        z_values.push_back(z);
        mean_z += z;
      }
      mean_z *= inv_n;

      double s2 = 0.0;
      for (double z : z_values) {
        const double dz = z - mean_z;
        s2 += dz * dz;
      }
      s2 /= static_cast<double>(n - 1);
      const double sigma_ratio =
        std::sqrt(s2 / (static_cast<double>(n) * result.denominator *
                         result.denominator));
      result.uncertainty[i] = std::abs(parameter_values_[i]) * sigma_ratio;
    }
  }

  return result;
}

ClutchSensitivityAccumulator::MethodResult
ClutchSensitivityAccumulator::compute_result() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return compute_result(batches_);
}

ClutchSensitivityAccumulator::MethodResult
ClutchSensitivityAccumulator::compute_cclutch_result() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return compute_cclutch_result(batches_);
}

vector<ClutchSensitivityAccumulator::BatchScore>
ClutchSensitivityAccumulator::collect_global_batches(
  int64_t& total_fission_sites, int64_t& total_scored_sites,
  int64_t& total_dropped_sites, int64_t& total_cclutch_events,
  int64_t& total_cclutch_scored_events,
  int64_t& total_cclutch_dropped_events,
  int64_t& total_cclutch_missing_source_events) const
{
  vector<BatchScore> local_batches;
  {
    std::lock_guard<std::mutex> lock(mutex_);
    local_batches = batches_;
    total_fission_sites = total_fission_sites_;
    total_scored_sites = total_scored_sites_;
    total_dropped_sites = total_dropped_sites_;
    total_cclutch_events = total_cclutch_events_;
    total_cclutch_scored_events = total_cclutch_scored_events_;
    total_cclutch_dropped_events = total_cclutch_dropped_events_;
    total_cclutch_missing_source_events =
      total_cclutch_missing_source_events_;
  }

#ifdef OPENMC_MPI
  int local_n = static_cast<int>(local_batches.size());
  int min_n = local_n;
  int max_n = local_n;
  MPI_Allreduce(&local_n, &min_n, 1, MPI_INT, MPI_MIN, mpi::intracomm);
  MPI_Allreduce(&local_n, &max_n, 1, MPI_INT, MPI_MAX, mpi::intracomm);
  if (min_n != max_n) {
    fatal_error("CLUTCH sensitivity batch counts differ across MPI ranks.");
  }

  const size_t n_params = n_parameters();
  const size_t block = 2 + 2 * n_params;
  vector<double> local(static_cast<size_t>(local_n) * block, 0.0);
  for (int b = 0; b < local_n; ++b) {
    size_t offset = static_cast<size_t>(b) * block;
    local[offset++] = local_batches[b].denominator;
    local[offset++] = local_batches[b].cclutch_denominator;
    for (double v : local_batches[b].numerator)
      local[offset++] = v;
    for (double v : local_batches[b].cclutch_numerator)
      local[offset++] = v;
  }

  vector<double> global;
  if (mpi::master) {
    global.assign(local.size(), 0.0);
  }
  if (!local.empty()) {
    MPI_Reduce(local.data(), mpi::master ? global.data() : nullptr,
      static_cast<int>(local.size()), MPI_DOUBLE, MPI_SUM, 0, mpi::intracomm);
  }

  int64_t local_counts[7] {total_fission_sites, total_scored_sites,
    total_dropped_sites, total_cclutch_events, total_cclutch_scored_events,
    total_cclutch_dropped_events, total_cclutch_missing_source_events};
  int64_t global_counts[7] {};
  MPI_Reduce(local_counts, global_counts, 7, MPI_INT64_T, MPI_SUM, 0,
    mpi::intracomm);

  if (!mpi::master) {
    return {};
  }
  total_fission_sites = global_counts[0];
  total_scored_sites = global_counts[1];
  total_dropped_sites = global_counts[2];
  total_cclutch_events = global_counts[3];
  total_cclutch_scored_events = global_counts[4];
  total_cclutch_dropped_events = global_counts[5];
  total_cclutch_missing_source_events = global_counts[6];

  vector<BatchScore> result = local_batches;
  for (auto& batch : result) {
    batch.numerator.assign(n_params, 0.0);
    batch.cclutch_numerator.assign(n_params, 0.0);
  }
  for (int b = 0; b < local_n; ++b) {
    size_t offset = static_cast<size_t>(b) * block;
    result[b].denominator = global[offset++];
    result[b].cclutch_denominator = global[offset++];
    for (int i = 0; i < n_params; ++i)
      result[b].numerator[i] = global[offset++];
    for (int i = 0; i < n_params; ++i)
      result[b].cclutch_numerator[i] = global[offset++];
  }
  return result;
#else
  return local_batches;
#endif
}

void ClutchSensitivityAccumulator::write_batch_matrix(hid_t group,
  const char* name, const vector<double>& flat, hsize_t n_batches) const
{
  hsize_t dims[] {n_batches, static_cast<hsize_t>(n_parameters())};
  write_dataset_lowlevel(group, 2, dims, name, H5TypeMap<double>::type_id,
    H5S_ALL, false, flat.data());
}

void ClutchSensitivityAccumulator::write_method_group(hid_t parent,
  const char* name, const MethodResult& result,
  const vector<BatchScore>& batches, bool use_cclutch) const
{
  hid_t group = create_group(parent, name);
  write_attribute(group, "available", static_cast<int>(result.available));
  write_attribute(group, "uses_fission_event_sites", use_cclutch ? 0 : 1);
  write_attribute(group, "uses_transfer_function", use_cclutch ? 1 : 0);
  write_attribute(group, "source_state_definition", "cell");
  if (use_cclutch) {
    write_attribute(group, "formula",
      "D=mean_b sum_source I*(source_cell) T_total(source_cell); "
      "N_x=mean_b sum_source I*(source_cell) T_x(source_cell)");
  } else {
    write_attribute(group, "formula",
      "D=mean_b sum_sites w_site I*(cell); "
      "N_x=mean_b sum_sites w_site I*(cell) (R_path_x+R_fission_x)");
  }
  write_dataset(group, "dlogk_dparameter", result.dlogk_dparameter);
  write_dataset(group, "sensitivity", result.sensitivity);
  write_dataset(group, "numerator", result.numerator);
  write_dataset(group, "denominator", vector<double> {result.denominator});
  write_dataset(group, "uncertainty", result.uncertainty);

  vector<int> batch_ids;
  vector<double> batch_denominator;
  vector<double> batch_numerator;
  batch_ids.reserve(batches.size());
  batch_denominator.reserve(batches.size());
  batch_numerator.reserve(batches.size() * n_parameters());
  for (const auto& batch : batches) {
    batch_ids.push_back(batch.batch_id);
    batch_denominator.push_back(
      use_cclutch ? batch.cclutch_denominator : batch.denominator);
    const auto& numerator =
      use_cclutch ? batch.cclutch_numerator : batch.numerator;
    batch_numerator.insert(
      batch_numerator.end(), numerator.begin(), numerator.end());
  }

  write_dataset(group, "batch_ids", batch_ids);
  write_dataset(group, "batch_denominator", batch_denominator);
  write_batch_matrix(
    group, "batch_numerator", batch_numerator, batches.size());

  H5Gclose(group);
}

void ClutchSensitivityAccumulator::write_to_file(
  const std::string& filename) const
{
  {
    std::lock_guard<std::mutex> lock(mutex_);
    if (batch_active_) {
      fatal_error("ClutchSensitivityAccumulator: cannot write while a batch is "
                  "active.");
    }
  }

  int64_t total_fission_sites = 0;
  int64_t total_scored_sites = 0;
  int64_t total_dropped_sites = 0;
  int64_t total_cclutch_events = 0;
  int64_t total_cclutch_scored_events = 0;
  int64_t total_cclutch_dropped_events = 0;
  int64_t total_cclutch_missing_source_events = 0;
  auto global_batches = collect_global_batches(
    total_fission_sites, total_scored_sites, total_dropped_sites,
    total_cclutch_events, total_cclutch_scored_events,
    total_cclutch_dropped_events, total_cclutch_missing_source_events);

  if (!mpi::master) {
    return;
  }

  auto result = compute_result(global_batches);
  auto cclutch_result = compute_cclutch_result(global_batches);
  if ((method_ == "hybrid" || method_ == "fclutch_fm") && !result.available) {
    fatal_error("CLUTCH sensitivity denominator is zero. Verify the spatial "
                "adjoint source and active fission-site scoring.");
  }
  if ((method_ == "hybrid" || method_ == "cclutch_history") &&
      !cclutch_result.available) {
    fatal_error("C-CLUTCH sensitivity denominator is zero. Verify active "
                "source birth recording and transfer-function scoring.");
  }

  hid_t file_id = file_open(filename, 'w');
  write_attribute(file_id, "filetype", "clutch_sensitivity");
  write_attribute(file_id, "version", "1.0");
  write_attribute(file_id, "sensitivity_type", "k_effective");
  write_attribute(file_id, "primary_method", "fclutch_fm");
  write_attribute(file_id, "method", method_);
  write_attribute(file_id, "adjoint_source", "fission_matrix.h5/adjoint_source");
  write_attribute(file_id, "cclutch_method",
    "cclutch_history is transfer-function C-CLUTCH in this version");

  hid_t params = create_group(file_id, "parameters");
  write_dataset(params, "ids", derivative_ids_);
  write_dataset(params, "variable", derivative_variables_);
  write_dataset(params, "material_id", derivative_material_ids_);
  write_dataset(params, "nuclide", derivative_nuclides_);
  write_dataset(params, "parameter_value", parameter_values_);
  H5Gclose(params);

  hid_t method_group = create_group(file_id, "method");
  if (method_ == "hybrid" || method_ == "fclutch_fm") {
    write_method_group(
      method_group, "fclutch_fm", result, global_batches, false);
  }
  if (method_ == "hybrid" || method_ == "cclutch_history") {
    write_method_group(
      method_group, "cclutch_history", cclutch_result, global_batches, true);
  }
  H5Gclose(method_group);

  hid_t diagnostics = create_group(file_id, "diagnostics");
  write_attribute(diagnostics, "total_fission_sites", total_fission_sites);
  write_attribute(diagnostics, "total_scored_sites", total_scored_sites);
  write_attribute(diagnostics, "total_dropped_sites", total_dropped_sites);
  write_attribute(diagnostics, "total_cclutch_events", total_cclutch_events);
  write_attribute(diagnostics, "total_cclutch_scored_events",
    total_cclutch_scored_events);
  write_attribute(diagnostics, "total_cclutch_dropped_events",
    total_cclutch_dropped_events);
  write_attribute(diagnostics, "total_cclutch_missing_source_events",
    total_cclutch_missing_source_events);
  write_attribute(diagnostics, "active_batches_scored",
    static_cast<int>(global_batches.size()));
  write_dataset(diagnostics, "grid_shape", grid_->shape());
  write_dataset(diagnostics, "grid_lower_left", grid_->origin());
  write_dataset(diagnostics, "grid_upper_right", grid_->upper_bound());
  std::array<double, 3> pitch {grid_->pitch(), grid_->pitch(), grid_->pitch()};
  write_dataset(diagnostics, "grid_pitch", pitch);
  H5Gclose(diagnostics);

  file_close(file_id);
}

int ClutchSensitivityAccumulator::position_to_index(const Position& r) const
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

} // namespace openmc
