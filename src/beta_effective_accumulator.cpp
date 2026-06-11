#include "openmc/beta_effective_accumulator.h"

#include "openmc/bank.h"
#include "openmc/clutch_ifp.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/mesh_init.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <unordered_set>

namespace openmc {

BetaEffectiveAccumulator::BetaEffectiveAccumulator(
  std::shared_ptr<SharedMeshGrid> grid, int n_materials)
  : grid_ {std::move(grid)}
{
  if (!grid_) {
    throw std::runtime_error(
      "BetaEffectiveAccumulator requires a valid SharedMeshGrid.");
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

void BetaEffectiveAccumulator::set_adjoint_source_spatial(
  const std::unordered_map<int64_t, double>& adjoint_source_spatial)
{
  bool has_positive = false;
  std::unordered_map<int64_t, double> cell_sum;
  std::unordered_map<int64_t, int> cell_count;
  for (const auto& [state, v] : adjoint_source_spatial) {
    if (state < 0 || state >= static_cast<int64_t>(n_source_states_)) {
      fatal_error("BetaEffectiveAccumulator: spatial adjoint source contains "
                  "an invalid cell-material state index.");
    }
    if (!std::isfinite(v)) {
      fatal_error("BetaEffectiveAccumulator: spatial adjoint source contains "
                  "non-finite values.");
    }
    has_positive = has_positive || v > 0.0;
    if (v > 0.0) {
      const int64_t cell = state / n_materials_;
      cell_sum[cell] += v;
      ++cell_count[cell];
    }
  }
  if (!has_positive) {
    fatal_error(
      "BetaEffectiveAccumulator: spatial adjoint source is all zero.");
  }

  std::unordered_map<int64_t, double> cell_fallback;
  cell_fallback.reserve(cell_sum.size());
  for (const auto& [cell, sum] : cell_sum) {
    const int count = cell_count[cell];
    if (count > 0) {
      cell_fallback[cell] = sum / static_cast<double>(count);
    }
  }

  std::lock_guard<std::mutex> lock(mutex_);
  adjoint_source_spatial_ = adjoint_source_spatial;
  adjoint_source_cell_fallback_ = std::move(cell_fallback);
  source_ready_ = true;
}

void BetaEffectiveAccumulator::begin_batch(int batch_id)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!source_ready_) {
    fatal_error("BetaEffectiveAccumulator: cannot begin an active beta_eff "
                "batch before the spatial adjoint source is loaded.");
  }
  if (batch_active_) {
    fatal_error("BetaEffectiveAccumulator: previous batch was not ended.");
  }
  current_batch_ = BatchScore {};
  current_batch_.batch_id = batch_id;
  current_source_states_.clear();
  current_source_counts_.clear();
  current_clutch_ifp_source_counts_.clear();
  current_fission_site_total_by_state_.clear();
  current_fission_site_delayed_by_state_.clear();
  current_clutch_ifp_response_by_state_.clear();
  current_cclutch_transfer_total_.clear();
  current_cclutch_transfer_lifetime_.clear();
  current_cclutch_transfer_emission_adjusted_lifetime_.clear();
  current_cclutch_transfer_delayed_.clear();
  batch_active_ = true;
}

void BetaEffectiveAccumulator::end_batch(int batch_id)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    return;
  }
  if (current_batch_.batch_id != batch_id) {
    fatal_error("BetaEffectiveAccumulator: end_batch called for a different "
                "batch than begin_batch.");
  }
  fold_current_cclutch_batch();
  fold_current_clutch_ifp_batch();
  batches_.push_back(current_batch_);
  batch_active_ = false;
  current_batch_ = BatchScore {};
  current_source_states_.clear();
  current_clutch_ifp_source_counts_.clear();
  current_fission_site_total_by_state_.clear();
  current_fission_site_delayed_by_state_.clear();
  current_clutch_ifp_response_by_state_.clear();
}

void BetaEffectiveAccumulator::record_source_birth(
  const Position& r, int64_t source_particle_id, int material_index,
  int64_t source_bank_index)
{
  if (source_particle_id < 0) {
    return;
  }

  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    return;
  }

  const int64_t state = source_state_index(r, material_index);
  if (state < 0) {
    ++total_invalid_source_states_;
    return;
  }

  auto [it, inserted] =
    current_source_states_.emplace(source_particle_id, state);
  if (inserted) {
    ++current_source_counts_[state];
  } else if (it->second != state) {
    --current_source_counts_[it->second];
    it->second = state;
    ++current_source_counts_[state];
  }

  if (source_bank_index < 0 ||
      source_bank_index >= static_cast<int64_t>(
                             simulation::clutch_ifp_source_state_bank.size())) {
    return;
  }

  const auto& states =
    simulation::clutch_ifp_source_state_bank[source_bank_index];
  if (states.size() != static_cast<size_t>(clutch_ifp_n_generation())) {
    return;
  }
  const int64_t ancestor_state = states.front();
  if (ancestor_state >= 0 &&
      ancestor_state < static_cast<int64_t>(n_source_states_)) {
    ++current_clutch_ifp_source_counts_[ancestor_state];
  }
}

void BetaEffectiveAccumulator::score_fission_site(const Position& r,
  double site_weight, int material_index, int delayed_group,
  double neutron_lifetime, double delayed_group_delay)
{
  if (!std::isfinite(site_weight) || site_weight <= 0.0) {
    return;
  }
  if (!std::isfinite(neutron_lifetime) || neutron_lifetime < 0.0) {
    fatal_error("BetaEffectiveAccumulator: fission-site neutron lifetime is "
                "negative or non-finite.");
  }
  if (!std::isfinite(delayed_group_delay) || delayed_group_delay < 0.0) {
    fatal_error("BetaEffectiveAccumulator: delayed-group mean delay is "
                "negative or non-finite.");
  }

  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    fatal_error("BetaEffectiveAccumulator: fission site scored outside an "
                "active beta_eff batch.");
  }
  if (!source_ready_) {
    fatal_error("BetaEffectiveAccumulator: adjoint source is not loaded.");
  }

  ++current_batch_.fission_sites;
  ++total_fission_sites_;

  const int64_t state = source_state_index(r, material_index);
  if (state < 0) {
    ++current_batch_.dropped_sites;
    ++total_dropped_sites_;
    ++total_invalid_fission_states_;
    return;
  }

  current_fission_site_total_by_state_[state] += site_weight;
  if (delayed_group > 0 && delayed_group <= N_DELAYED_GROUPS) {
    current_fission_site_delayed_by_state_[state][delayed_group - 1] +=
      site_weight;
  }

  auto importance_it = adjoint_source_spatial_.find(state);
  double importance = 0.0;
  bool used_fallback = false;
  if (importance_it == adjoint_source_spatial_.end() ||
      importance_it->second <= 0.0) {
    const int64_t cell = state / n_materials_;
    auto fallback_it = adjoint_source_cell_fallback_.find(cell);
    if (fallback_it == adjoint_source_cell_fallback_.end() ||
        fallback_it->second <= 0.0) {
      ++current_batch_.dropped_sites;
      ++total_dropped_sites_;
      return;
    }
    importance = fallback_it->second;
    used_fallback = true;
  } else {
    importance = importance_it->second;
  }

  const double score_weight = site_weight * importance;
  current_batch_.denominator += score_weight;
  current_batch_.lifetime_numerator += score_weight * neutron_lifetime;
  current_batch_.emission_adjusted_lifetime_numerator +=
    score_weight * neutron_lifetime;

  if (delayed_group < 0 || delayed_group > N_DELAYED_GROUPS) {
    ++current_batch_.invalid_delayed_group_sites;
    ++total_invalid_delayed_group_sites_;
  } else if (delayed_group > 0) {
    current_batch_.numerator[delayed_group - 1] += score_weight;
    current_batch_.emission_adjusted_lifetime_numerator +=
      score_weight * delayed_group_delay;
  }

  ++current_batch_.scored_sites;
  ++total_scored_sites_;
  if (used_fallback) {
    ++current_batch_.fallback_sites;
    ++total_fallback_sites_;
  }
}

void BetaEffectiveAccumulator::score_cclutch_fission_event(const Position& r,
  int64_t source_particle_id, double total_contribution,
  const std::array<double, N_DELAYED_GROUPS>& delayed_contributions,
  double neutron_lifetime,
  const std::array<double, N_DELAYED_GROUPS>& delayed_group_delays)
{
  if (!std::isfinite(total_contribution) || total_contribution <= 0.0) {
    return;
  }
  if (!std::isfinite(neutron_lifetime) || neutron_lifetime < 0.0) {
    fatal_error("BetaEffectiveAccumulator: C-CLUTCH neutron lifetime is "
                "negative or non-finite.");
  }

  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    fatal_error("BetaEffectiveAccumulator: C-CLUTCH event scored outside an "
                "active beta_eff batch.");
  }
  if (!source_ready_) {
    fatal_error("BetaEffectiveAccumulator: adjoint source is not loaded.");
  }

  ++current_batch_.cclutch_events;
  ++total_cclutch_events_;

  if (position_to_index(r) < 0) {
    ++current_batch_.cclutch_dropped_events;
    ++total_cclutch_dropped_events_;
    return;
  }

  auto it = current_source_states_.find(source_particle_id);
  if (it == current_source_states_.end()) {
    ++current_batch_.cclutch_missing_source_events;
    ++total_cclutch_missing_source_events_;
    return;
  }

  const int64_t source_state = it->second;
  if (source_state < 0 ||
      source_state >= static_cast<int64_t>(n_source_states_)) {
    ++current_batch_.cclutch_dropped_events;
    ++total_cclutch_dropped_events_;
    return;
  }

  current_cclutch_transfer_total_[source_state] += total_contribution;
  current_cclutch_transfer_lifetime_[source_state] +=
    total_contribution * neutron_lifetime;
  current_cclutch_transfer_emission_adjusted_lifetime_[source_state] +=
    total_contribution * neutron_lifetime;
  for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
    const double delayed = delayed_contributions[k];
    if (std::isfinite(delayed) && delayed > 0.0) {
      const double delay = delayed_group_delays[k];
      if (!std::isfinite(delay) || delay <= 0.0) {
        fatal_error("BetaEffectiveAccumulator: C-CLUTCH delayed response has "
                    "no positive delayed-group mean delay.");
      }
      current_cclutch_transfer_delayed_[source_state][k] += delayed;
      current_cclutch_transfer_emission_adjusted_lifetime_[source_state] +=
        delayed * delay;
    }
  }

  ++current_batch_.cclutch_scored_events;
  ++total_cclutch_scored_events_;
}

void BetaEffectiveAccumulator::score_ifp_ancestry_event(
  double fission_weight, int64_t source_bank_index)
{
  if (!std::isfinite(fission_weight) || fission_weight <= 0.0) {
    return;
  }

  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    return;
  }

  ++current_batch_.ifp_ancestry_events;
  ++total_ifp_ancestry_events_;

  if (source_bank_index < 0 ||
      source_bank_index >= static_cast<int64_t>(
                             simulation::clutch_ifp_source_state_bank.size()) ||
      source_bank_index >= static_cast<int64_t>(
                             simulation::clutch_ifp_source_delayed_group_bank
                               .size())) {
    ++current_batch_.ifp_ancestry_incomplete_events;
    ++total_ifp_ancestry_incomplete_events_;
    return;
  }

  const auto& states =
    simulation::clutch_ifp_source_state_bank[source_bank_index];
  const auto& delayed_groups =
    simulation::clutch_ifp_source_delayed_group_bank[source_bank_index];
  const int n_generation = clutch_ifp_n_generation();
  if (states.size() != static_cast<size_t>(n_generation) ||
      delayed_groups.size() != static_cast<size_t>(n_generation)) {
    ++current_batch_.ifp_ancestry_incomplete_events;
    ++total_ifp_ancestry_incomplete_events_;
    return;
  }

  const int64_t ancestor_state = states.front();
  if (ancestor_state < 0 ||
      ancestor_state >= static_cast<int64_t>(n_source_states_)) {
    ++current_batch_.ifp_ancestry_incomplete_events;
    ++total_ifp_ancestry_incomplete_events_;
    return;
  }

  current_clutch_ifp_response_by_state_[ancestor_state] += fission_weight;

  const int delayed_group = delayed_groups.front();
  current_batch_.ifp_ancestry_denominator += fission_weight;
  if (delayed_group < 0 || delayed_group > N_DELAYED_GROUPS) {
    ++current_batch_.ifp_ancestry_invalid_delayed_group_events;
    ++total_ifp_ancestry_invalid_delayed_group_events_;
  } else if (delayed_group > 0) {
    current_batch_.ifp_ancestry_numerator[delayed_group - 1] +=
      fission_weight;
  }

  ++current_batch_.ifp_ancestry_scored_events;
  ++total_ifp_ancestry_scored_events_;
}

void BetaEffectiveAccumulator::fold_current_cclutch_batch()
{
  std::unordered_set<int64_t> source_states;
  source_states.reserve(current_cclutch_transfer_total_.size() +
                        current_cclutch_transfer_delayed_.size());
  for (const auto& [source_state, total] : current_cclutch_transfer_total_) {
    source_states.insert(source_state);
  }
  for (const auto& [source_state, delayed] : current_cclutch_transfer_delayed_) {
    source_states.insert(source_state);
  }

  for (int64_t source_state : source_states) {
    auto total_it = current_cclutch_transfer_total_.find(source_state);
    const double total =
      total_it != current_cclutch_transfer_total_.end() ? total_it->second :
                                                          0.0;
    auto delayed_it = current_cclutch_transfer_delayed_.find(source_state);
    bool has_delayed = false;
    if (delayed_it != current_cclutch_transfer_delayed_.end()) {
      for (double delayed : delayed_it->second) {
        has_delayed = has_delayed || delayed > 0.0;
      }
    }
    if (total <= 0.0 && !has_delayed) {
      continue;
    }

    auto count_it = current_source_counts_.find(source_state);
    const int source_count =
      count_it != current_source_counts_.end() ? count_it->second : 0;
    if (source_count <= 0) {
      fatal_error("BetaEffectiveAccumulator: C-CLUTCH transfer response exists "
                  "for a source state with zero source count.");
    }

    auto importance_it = adjoint_source_spatial_.find(source_state);
    if (importance_it == adjoint_source_spatial_.end() ||
        importance_it->second <= 0.0) {
      continue;
    }

    const double source_weight =
      importance_it->second / static_cast<double>(source_count);
    current_batch_.cclutch_denominator += source_weight * total;
    auto lifetime_it = current_cclutch_transfer_lifetime_.find(source_state);
    current_batch_.cclutch_lifetime_numerator +=
      source_weight *
      (lifetime_it != current_cclutch_transfer_lifetime_.end() ?
          lifetime_it->second :
          0.0);
    auto adjusted_it =
      current_cclutch_transfer_emission_adjusted_lifetime_.find(source_state);
    current_batch_.cclutch_emission_adjusted_lifetime_numerator +=
      source_weight *
      (adjusted_it !=
          current_cclutch_transfer_emission_adjusted_lifetime_.end() ?
          adjusted_it->second :
          0.0);
    if (delayed_it != current_cclutch_transfer_delayed_.end()) {
      for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
        current_batch_.cclutch_numerator[k] +=
          source_weight * delayed_it->second[k];
      }
    }
  }
}

void BetaEffectiveAccumulator::fold_current_clutch_ifp_batch()
{
  if (current_fission_site_total_by_state_.empty() ||
      current_clutch_ifp_response_by_state_.empty()) {
    return;
  }

  std::unordered_map<int64_t, double> ifp_importance;
  ifp_importance.reserve(current_clutch_ifp_response_by_state_.size());
  for (const auto& [state, response] : current_clutch_ifp_response_by_state_) {
    if (!std::isfinite(response) || response <= 0.0) {
      continue;
    }
    auto count_it = current_clutch_ifp_source_counts_.find(state);
    const int source_count =
      count_it != current_clutch_ifp_source_counts_.end() ? count_it->second :
                                                            0;
    if (source_count > 0) {
      ifp_importance[state] = response / static_cast<double>(source_count);
    }
  }

  if (ifp_importance.empty()) {
    return;
  }

  for (const auto& [state, total] : current_fission_site_total_by_state_) {
    if (!std::isfinite(total) || total <= 0.0) {
      continue;
    }
    auto importance_it = ifp_importance.find(state);
    if (importance_it == ifp_importance.end() ||
        importance_it->second <= 0.0) {
      continue;
    }

    const double weight = total * importance_it->second;
    current_batch_.clutch_ifp_denominator += weight;

    auto delayed_it = current_fission_site_delayed_by_state_.find(state);
    if (delayed_it != current_fission_site_delayed_by_state_.end()) {
      for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
        current_batch_.clutch_ifp_numerator[k] +=
          delayed_it->second[k] * importance_it->second;
      }
    }
  }
}

BetaEffectiveAccumulator::MethodResult
BetaEffectiveAccumulator::compute_result() const
{
  return compute_result(MethodKind::FClutch);
}

BetaEffectiveAccumulator::MethodResult
BetaEffectiveAccumulator::compute_cclutch_result() const
{
  return compute_result(MethodKind::CClutch);
}

BetaEffectiveAccumulator::MethodResult
BetaEffectiveAccumulator::compute_ifp_ancestry_result() const
{
  return compute_result(MethodKind::IfpAncestry);
}

BetaEffectiveAccumulator::MethodResult
BetaEffectiveAccumulator::compute_clutch_ifp_result() const
{
  return compute_result(MethodKind::ClutchIfp);
}

BetaEffectiveAccumulator::GenerationTimeResult
BetaEffectiveAccumulator::compute_generation_time_result() const
{
  return compute_generation_time_result(false);
}

BetaEffectiveAccumulator::GenerationTimeResult
BetaEffectiveAccumulator::compute_cclutch_generation_time_result() const
{
  return compute_generation_time_result(true);
}

BetaEffectiveAccumulator::MethodResult
BetaEffectiveAccumulator::compute_result(MethodKind method) const
{
  std::lock_guard<std::mutex> lock(mutex_);

  MethodResult result;
  const size_t n = batches_.size();
  if (n == 0) {
    return result;
  }

  double sum_denominator = 0.0;
  std::array<double, N_DELAYED_GROUPS> sum_numerator {};
  for (const auto& batch : batches_) {
    switch (method) {
    case MethodKind::FClutch:
      sum_denominator += batch.denominator;
      break;
    case MethodKind::CClutch:
      sum_denominator += batch.cclutch_denominator;
      break;
    case MethodKind::IfpAncestry:
      sum_denominator += batch.ifp_ancestry_denominator;
      break;
    case MethodKind::ClutchIfp:
      sum_denominator += batch.clutch_ifp_denominator;
      break;
    }
    for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
      switch (method) {
      case MethodKind::FClutch:
        sum_numerator[k] += batch.numerator[k];
        break;
      case MethodKind::CClutch:
        sum_numerator[k] += batch.cclutch_numerator[k];
        break;
      case MethodKind::IfpAncestry:
        sum_numerator[k] += batch.ifp_ancestry_numerator[k];
        break;
      case MethodKind::ClutchIfp:
        sum_numerator[k] += batch.clutch_ifp_numerator[k];
        break;
      }
    }
  }

  if (sum_denominator <= 0.0) {
    return result;
  }

  const double inv_n = 1.0 / static_cast<double>(n);
  const double mean_denominator = sum_denominator * inv_n;
  result.denominator = mean_denominator;
  result.available = true;

  double total_numerator = 0.0;
  for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
    const double mean_numerator = sum_numerator[k] * inv_n;
    result.numerators[k] = mean_numerator;
    result.beta_i[k] = mean_numerator / mean_denominator;
    result.beta_total += result.beta_i[k];
    total_numerator += mean_numerator;
  }

  if (n > 1) {
    for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
      const double beta = result.beta_i[k];
      double mean_z = 0.0;
      std::vector<double> z_values;
      z_values.reserve(n);
      for (const auto& batch : batches_) {
        double numerator = 0.0;
        double denominator = 0.0;
        switch (method) {
        case MethodKind::FClutch:
          numerator = batch.numerator[k];
          denominator = batch.denominator;
          break;
        case MethodKind::CClutch:
          numerator = batch.cclutch_numerator[k];
          denominator = batch.cclutch_denominator;
          break;
        case MethodKind::IfpAncestry:
          numerator = batch.ifp_ancestry_numerator[k];
          denominator = batch.ifp_ancestry_denominator;
          break;
        case MethodKind::ClutchIfp:
          numerator = batch.clutch_ifp_numerator[k];
          denominator = batch.clutch_ifp_denominator;
          break;
        }
        const double z = numerator - beta * denominator;
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
      result.uncertainty[k] =
        std::sqrt(s2 / (static_cast<double>(n) * mean_denominator *
                         mean_denominator));
    }

    const double beta_total =
      mean_denominator > 0.0 ? total_numerator / mean_denominator : 0.0;
    double mean_z = 0.0;
    std::vector<double> z_values;
    z_values.reserve(n);
    for (const auto& batch : batches_) {
      const auto& numerator =
        method == MethodKind::CClutch ? batch.cclutch_numerator :
        method == MethodKind::IfpAncestry ? batch.ifp_ancestry_numerator :
        method == MethodKind::ClutchIfp ? batch.clutch_ifp_numerator :
                                            batch.numerator;
      const double numerator_sum =
        std::accumulate(numerator.begin(), numerator.end(), 0.0);
      const double denominator =
        method == MethodKind::CClutch ? batch.cclutch_denominator :
        method == MethodKind::IfpAncestry ? batch.ifp_ancestry_denominator :
        method == MethodKind::ClutchIfp ? batch.clutch_ifp_denominator :
                                            batch.denominator;
      const double z = numerator_sum - beta_total * denominator;
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
    result.beta_total_uncertainty =
      std::sqrt(s2 / (static_cast<double>(n) * mean_denominator *
                       mean_denominator));
  }

  return result;
}

BetaEffectiveAccumulator::GenerationTimeResult
BetaEffectiveAccumulator::compute_generation_time_result(
  bool use_cclutch) const
{
  std::lock_guard<std::mutex> lock(mutex_);

  GenerationTimeResult result;
  const size_t n = batches_.size();
  if (n == 0) {
    return result;
  }

  double sum_denominator = 0.0;
  double sum_lifetime = 0.0;
  double sum_emission_adjusted = 0.0;
  for (const auto& batch : batches_) {
    sum_denominator +=
      use_cclutch ? batch.cclutch_denominator : batch.denominator;
    sum_lifetime += use_cclutch ? batch.cclutch_lifetime_numerator
                                : batch.lifetime_numerator;
    sum_emission_adjusted +=
      use_cclutch ? batch.cclutch_emission_adjusted_lifetime_numerator
                  : batch.emission_adjusted_lifetime_numerator;
  }

  if (sum_denominator <= 0.0) {
    return result;
  }

  const double inv_n = 1.0 / static_cast<double>(n);
  const double mean_denominator = sum_denominator * inv_n;
  const double mean_lifetime = sum_lifetime * inv_n;
  const double mean_emission_adjusted = sum_emission_adjusted * inv_n;
  result.available = true;
  result.denominator = mean_denominator;
  result.lifetime_numerator = mean_lifetime;
  result.emission_adjusted_lifetime_numerator = mean_emission_adjusted;
  result.transport_lifetime = mean_lifetime / mean_denominator;
  result.emission_adjusted_lifetime =
    mean_emission_adjusted / mean_denominator;

  if (n > 1) {
    auto ratio_uncertainty = [&](double ratio, bool adjusted) {
      double mean_z = 0.0;
      std::vector<double> z_values;
      z_values.reserve(n);
      for (const auto& batch : batches_) {
        const double numerator =
          adjusted ? (use_cclutch
                         ? batch.cclutch_emission_adjusted_lifetime_numerator
                         : batch.emission_adjusted_lifetime_numerator)
                   : (use_cclutch ? batch.cclutch_lifetime_numerator
                                  : batch.lifetime_numerator);
        const double denominator =
          use_cclutch ? batch.cclutch_denominator : batch.denominator;
        const double z = numerator - ratio * denominator;
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
      return std::sqrt(s2 / (static_cast<double>(n) * mean_denominator *
                              mean_denominator));
    };

    result.transport_lifetime_uncertainty =
      ratio_uncertainty(result.transport_lifetime, false);
    result.emission_adjusted_lifetime_uncertainty =
      ratio_uncertainty(result.emission_adjusted_lifetime, true);
  }

  return result;
}

size_t BetaEffectiveAccumulator::n_batches() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return batches_.size();
}

size_t BetaEffectiveAccumulator::scored_sites() const
{
  std::lock_guard<std::mutex> lock(mutex_);
  return static_cast<size_t>(total_scored_sites_);
}

int BetaEffectiveAccumulator::position_to_index(const Position& r) const
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

int64_t BetaEffectiveAccumulator::source_state_index(
  const Position& r, int material_index) const
{
  if (material_index < 0 || material_index >= n_materials_)
    return -1;

  const int cell = position_to_index(r);
  if (cell < 0)
    return -1;

  return static_cast<int64_t>(cell) * n_materials_ + material_index;
}

void BetaEffectiveAccumulator::write_batch_matrix(hid_t group,
  const char* name, const std::vector<double>& flat, hsize_t n_batches) const
{
  hsize_t dims[] {n_batches, static_cast<hsize_t>(N_DELAYED_GROUPS)};
  write_dataset_lowlevel(group, 2, dims, name, H5TypeMap<double>::type_id,
    H5S_ALL, false, flat.data());
}

void BetaEffectiveAccumulator::write_method_group(hid_t parent,
  const char* name, const MethodResult& result, bool use_cclutch) const
{
  hid_t group = create_group(parent, name);
  write_attribute(group, "available", static_cast<int>(result.available));
  write_attribute(group, "theory_reference",
    use_cclutch ? "Qiu2016_C_CLUTCH_transfer_function"
                : "Qiu2016_F_CLUTCH_Eq31_Eq43");
  write_attribute(group, "uses_fission_event_sites", use_cclutch ? 0 : 1);
  write_attribute(group, "uses_transfer_function", use_cclutch ? 1 : 0);
  write_attribute(group, "uses_birth_energy_importance", 0);
  write_attribute(group, "source_state_definition", "cell_material");
  if (!use_cclutch) {
    write_attribute(group, "missing_importance_fallback",
      "cell_material_mean_importance_in_same_cell");
  }
  if (use_cclutch) {
    write_attribute(group, "formula",
      "D=mean_b sum_source I*(source_cell,source_material) "
      "T_total(source_cell,source_material); "
      "N_k=mean_b sum_source I*(source_cell,source_material) "
      "T_delayed_k(source_cell,source_material)");
  } else {
    write_attribute(group, "formula",
      "D=mean_b sum_sites w_site I*(cell,material); "
      "N_k=mean_b sum_delayed_k w_site I*(cell,material)");
  }

  write_dataset(
    group, "beta_i", std::vector<double>(result.beta_i.begin(), result.beta_i.end()));
  write_dataset(group, "beta_total", std::vector<double> {result.beta_total});
  write_dataset(group, "numerator",
    std::vector<double>(result.numerators.begin(), result.numerators.end()));
  write_dataset(group, "denominator", std::vector<double> {result.denominator});
  write_dataset(group, "uncertainty",
    std::vector<double>(result.uncertainty.begin(), result.uncertainty.end()));
  write_dataset(group, "beta_total_uncertainty",
    std::vector<double> {result.beta_total_uncertainty});

  std::vector<int> batch_ids;
  std::vector<double> batch_denominator;
  std::vector<double> batch_numerator;
  batch_ids.reserve(batches_.size());
  batch_denominator.reserve(batches_.size());
  batch_numerator.reserve(
    batches_.size() * static_cast<size_t>(N_DELAYED_GROUPS));

  for (const auto& batch : batches_) {
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
    group, "batch_numerator", batch_numerator, batches_.size());

  H5Gclose(group);
}

void BetaEffectiveAccumulator::write_ifp_ancestry_method_group(
  hid_t parent, const MethodResult& result) const
{
  hid_t group = create_group(parent, "ifp_ancestry");
  write_attribute(group, "available", static_cast<int>(result.available));
  write_attribute(group, "theory_reference", "Hurwitz1964_IFP_ancestry");
  write_attribute(group, "uses_ifp", 1);
  write_attribute(group, "ifp_n_generation", clutch_ifp_n_generation());
  write_attribute(group, "uses_fission_event_sites", 0);
  write_attribute(group, "uses_transfer_function", 0);
  write_attribute(group, "uses_birth_energy_importance", 0);
  write_attribute(group, "source_state_definition", "cell_material");
  write_attribute(group, "formula",
    "D=mean_b sum_fission_events w_event; "
    "N_k=mean_b sum_fission_events w_event "
    "1[ancestor_delayed_group(k,N_gen)]");

  write_dataset(
    group, "beta_i", std::vector<double>(result.beta_i.begin(), result.beta_i.end()));
  write_dataset(group, "beta_total", std::vector<double> {result.beta_total});
  write_dataset(group, "numerator",
    std::vector<double>(result.numerators.begin(), result.numerators.end()));
  write_dataset(group, "denominator", std::vector<double> {result.denominator});
  write_dataset(group, "uncertainty",
    std::vector<double>(result.uncertainty.begin(), result.uncertainty.end()));
  write_dataset(group, "beta_total_uncertainty",
    std::vector<double> {result.beta_total_uncertainty});

  std::vector<int> batch_ids;
  std::vector<double> batch_denominator;
  std::vector<double> batch_numerator;
  batch_ids.reserve(batches_.size());
  batch_denominator.reserve(batches_.size());
  batch_numerator.reserve(
    batches_.size() * static_cast<size_t>(N_DELAYED_GROUPS));

  for (const auto& batch : batches_) {
    batch_ids.push_back(batch.batch_id);
    batch_denominator.push_back(batch.ifp_ancestry_denominator);
    batch_numerator.insert(batch_numerator.end(),
      batch.ifp_ancestry_numerator.begin(),
      batch.ifp_ancestry_numerator.end());
  }

  write_dataset(group, "batch_ids", batch_ids);
  write_dataset(group, "batch_denominator", batch_denominator);
  write_batch_matrix(
    group, "batch_numerator", batch_numerator, batches_.size());

  H5Gclose(group);
}

void BetaEffectiveAccumulator::write_clutch_ifp_method_group(
  hid_t parent, const MethodResult& result) const
{
  hid_t group = create_group(parent, "clutch_ifp");
  write_attribute(group, "available", static_cast<int>(result.available));
  write_attribute(group, "theory_reference",
    "CLUTCH_with_IFP_derived_cell_material_importance");
  write_attribute(group, "uses_ifp", 1);
  write_attribute(group, "ifp_n_generation", clutch_ifp_n_generation());
  write_attribute(group, "uses_fission_event_sites", 1);
  write_attribute(group, "uses_transfer_function", 0);
  write_attribute(group, "uses_birth_energy_importance", 0);
  write_attribute(group, "source_state_definition", "cell_material");
  write_attribute(group, "importance_definition",
    "I_ifp(state)=sum_descendant_fission_event_weight(state)/"
    "source_count_with_ancestor_state(state)");
  write_attribute(group, "formula",
    "D=mean_b sum_sites w_site I_ifp(cell,material); "
    "N_k=mean_b sum_delayed_k w_site I_ifp(cell,material)");

  write_dataset(
    group, "beta_i", std::vector<double>(result.beta_i.begin(), result.beta_i.end()));
  write_dataset(group, "beta_total", std::vector<double> {result.beta_total});
  write_dataset(group, "numerator",
    std::vector<double>(result.numerators.begin(), result.numerators.end()));
  write_dataset(group, "denominator", std::vector<double> {result.denominator});
  write_dataset(group, "uncertainty",
    std::vector<double>(result.uncertainty.begin(), result.uncertainty.end()));
  write_dataset(group, "beta_total_uncertainty",
    std::vector<double> {result.beta_total_uncertainty});

  std::vector<int> batch_ids;
  std::vector<double> batch_denominator;
  std::vector<double> batch_numerator;
  batch_ids.reserve(batches_.size());
  batch_denominator.reserve(batches_.size());
  batch_numerator.reserve(
    batches_.size() * static_cast<size_t>(N_DELAYED_GROUPS));

  for (const auto& batch : batches_) {
    batch_ids.push_back(batch.batch_id);
    batch_denominator.push_back(batch.clutch_ifp_denominator);
    batch_numerator.insert(batch_numerator.end(),
      batch.clutch_ifp_numerator.begin(), batch.clutch_ifp_numerator.end());
  }

  write_dataset(group, "batch_ids", batch_ids);
  write_dataset(group, "batch_denominator", batch_denominator);
  write_batch_matrix(
    group, "batch_numerator", batch_numerator, batches_.size());

  H5Gclose(group);
}

void BetaEffectiveAccumulator::write_generation_time_method_group(hid_t parent,
  const char* name, const GenerationTimeResult& result, bool use_cclutch) const
{
  hid_t group = create_group(parent, name);
  write_attribute(group, "available", static_cast<int>(result.available));
  write_attribute(group, "uses_ifp", 0);
  write_attribute(group, "units", "s");
  write_attribute(group, "source_state_definition", "cell_material");
  write_attribute(group, "uses_transfer_function", use_cclutch ? 1 : 0);
  write_attribute(group, "delay_treatment", "group_mean_1_over_lambda");
  write_attribute(group, "transport_lifetime_definition",
    "Particle::lifetime() from source birth to fission event");
  write_attribute(group, "emission_adjusted_lifetime_definition",
    "transport_lifetime plus delayed group mean precursor delay for delayed "
    "response only");
  if (use_cclutch) {
    write_attribute(group, "formula",
      "D=mean_b sum_source I*(source_cell,source_material) "
      "T_total(source_cell,source_material); "
      "L=mean_b sum_source I*(source_cell,source_material) "
      "T_lifetime(source_cell,source_material)");
  } else {
    write_attribute(group, "formula",
      "D=mean_b sum_sites w_site I*(cell,material); "
      "L=mean_b sum_sites w_site I*(cell,material) lifetime");
  }

  write_dataset(
    group, "transport_lifetime", std::vector<double> {result.transport_lifetime});
  write_dataset(group, "emission_adjusted_lifetime",
    std::vector<double> {result.emission_adjusted_lifetime});
  write_dataset(group, "denominator", std::vector<double> {result.denominator});
  write_dataset(group, "lifetime_numerator",
    std::vector<double> {result.lifetime_numerator});
  write_dataset(group, "emission_adjusted_lifetime_numerator",
    std::vector<double> {result.emission_adjusted_lifetime_numerator});
  write_dataset(group, "uncertainty",
    std::vector<double> {result.transport_lifetime_uncertainty});
  write_dataset(group, "emission_adjusted_uncertainty",
    std::vector<double> {result.emission_adjusted_lifetime_uncertainty});

  std::vector<int> batch_ids;
  std::vector<double> batch_denominator;
  std::vector<double> batch_lifetime_numerator;
  std::vector<double> batch_emission_adjusted_lifetime_numerator;
  batch_ids.reserve(batches_.size());
  batch_denominator.reserve(batches_.size());
  batch_lifetime_numerator.reserve(batches_.size());
  batch_emission_adjusted_lifetime_numerator.reserve(batches_.size());

  for (const auto& batch : batches_) {
    batch_ids.push_back(batch.batch_id);
    batch_denominator.push_back(
      use_cclutch ? batch.cclutch_denominator : batch.denominator);
    batch_lifetime_numerator.push_back(
      use_cclutch ? batch.cclutch_lifetime_numerator
                  : batch.lifetime_numerator);
    batch_emission_adjusted_lifetime_numerator.push_back(
      use_cclutch ? batch.cclutch_emission_adjusted_lifetime_numerator
                  : batch.emission_adjusted_lifetime_numerator);
  }

  write_dataset(group, "batch_ids", batch_ids);
  write_dataset(group, "batch_denominator", batch_denominator);
  write_dataset(group, "batch_lifetime_numerator", batch_lifetime_numerator);
  write_dataset(group, "batch_emission_adjusted_lifetime_numerator",
    batch_emission_adjusted_lifetime_numerator);

  H5Gclose(group);
}

void BetaEffectiveAccumulator::write_to_file(const std::string& filename) const
{
  if (batch_active_) {
    fatal_error(
      "BetaEffectiveAccumulator: cannot write beta_eff while a batch is active.");
  }

  auto result = compute_result();
  auto cclutch_result = compute_cclutch_result();
  auto ifp_ancestry_result = compute_ifp_ancestry_result();
  auto clutch_ifp_result = compute_clutch_ifp_result();
  if (!result.available) {
    fatal_error("F-CLUTCH spatial beta_eff denominator is zero. Verify the "
                "spatial adjoint source and active fission-site scoring.");
  }
  if (!cclutch_result.available) {
    fatal_error("C-CLUTCH beta_eff denominator is zero. Verify active source "
                "birth recording and C-CLUTCH fission event scoring.");
  }
  if (total_invalid_delayed_group_sites_ > 0) {
    warning("F-CLUTCH beta_eff ignored numerator contributions for " +
            std::to_string(total_invalid_delayed_group_sites_) +
            " fission sites with delayed_group outside [0, 8].");
  }
  if (total_cclutch_missing_source_events_ > 0) {
    warning("C-CLUTCH beta_eff ignored " +
            std::to_string(total_cclutch_missing_source_events_) +
            " fission response events without a recorded source cell.");
  }
  if (total_ifp_ancestry_incomplete_events_ > 0) {
    warning("IFP-ancestry beta_eff ignored " +
            std::to_string(total_ifp_ancestry_incomplete_events_) +
            " active fission events without a complete ancestry chain.");
  }

  hid_t file_id = file_open(filename, 'w');
  write_attribute(file_id, "filetype", "beta_effective");
  write_attribute(file_id, "version", "4.1");
  write_attribute(
    file_id, "description", "Effective delayed neutron fraction beta_eff");

  write_dataset(file_id, "beta_i",
    std::vector<double>(result.beta_i.begin(), result.beta_i.end()));
  write_dataset(file_id, "beta_total", std::vector<double> {result.beta_total});
  write_dataset(file_id, "uncertainty",
    std::vector<double>(result.uncertainty.begin(), result.uncertainty.end()));

  hid_t metadata = create_group(file_id, "metadata");
  write_attribute(metadata, "primary_method", "fclutch_cell_material");
  write_attribute(metadata, "theory_reference", "Qiu2016_F_CLUTCH_Eq31_Eq43");
  write_attribute(metadata, "uses_fission_event_sites", 1);
  write_attribute(metadata, "has_cclutch_method", 1);
  write_attribute(metadata, "has_ifp_ancestry_method", 1);
  write_attribute(metadata, "has_clutch_ifp_method", 1);
  write_attribute(metadata, "ifp_ancestry_n_generation",
    clutch_ifp_n_generation());
  write_attribute(metadata, "fclutch_missing_importance_fallback",
    "cell_material_mean_importance_in_same_cell");
  write_attribute(metadata, "uses_birth_energy_importance", 0);
  write_attribute(metadata, "source_state_definition", "cell_material");
  write_attribute(metadata, "state_indexing",
    "state=cell*n_materials+material_index");
  write_attribute(metadata, "n_delayed_groups", N_DELAYED_GROUPS);
  write_attribute(metadata, "n_materials", n_materials_);
  write_attribute(
    metadata, "n_source_states", static_cast<int64_t>(n_source_states_));
  write_attribute(metadata, "batchwise_ratio_uncertainty", 1);
  write_attribute(metadata, "adjoint_source",
    "fission_matrix.h5/adjoint_source generated by I*=(1/k)F^T I* on the "
    "energy-integrated cell-material fission matrix");
  write_dataset(metadata, "grid_shape", grid_->shape());
  write_dataset(metadata, "grid_lower_left", grid_->origin());
  write_dataset(metadata, "grid_upper_right", grid_->upper_bound());
  std::array<double, 3> pitch {grid_->pitch(), grid_->pitch(), grid_->pitch()};
  write_dataset(metadata, "grid_pitch", pitch);
  std::vector<int> material_ids;
  material_ids.reserve(n_materials_);
  for (int m = 0; m < n_materials_; ++m) {
    if (m < static_cast<int>(model::materials.size()) && model::materials[m]) {
      material_ids.push_back(model::materials[m]->id());
    } else {
      material_ids.push_back(m);
    }
  }
  write_dataset(metadata, "material_ids", material_ids);
  H5Gclose(metadata);

  hid_t diagnostics = create_group(file_id, "diagnostics");
  write_attribute(diagnostics, "total_fission_sites", total_fission_sites_);
  write_attribute(diagnostics, "total_scored_sites", total_scored_sites_);
  write_attribute(diagnostics, "total_dropped_sites", total_dropped_sites_);
  write_attribute(diagnostics, "total_fallback_sites", total_fallback_sites_);
  write_attribute(diagnostics, "total_invalid_delayed_group_sites",
    total_invalid_delayed_group_sites_);
  write_attribute(diagnostics, "total_invalid_source_states",
    total_invalid_source_states_);
  write_attribute(diagnostics, "total_invalid_fission_states",
    total_invalid_fission_states_);
  write_attribute(diagnostics, "total_cclutch_events", total_cclutch_events_);
  write_attribute(diagnostics, "total_cclutch_scored_events",
    total_cclutch_scored_events_);
  write_attribute(diagnostics, "total_cclutch_dropped_events",
    total_cclutch_dropped_events_);
  write_attribute(diagnostics, "total_cclutch_missing_source_events",
    total_cclutch_missing_source_events_);
  write_attribute(diagnostics, "total_ifp_ancestry_events",
    total_ifp_ancestry_events_);
  write_attribute(diagnostics, "total_ifp_ancestry_scored_events",
    total_ifp_ancestry_scored_events_);
  write_attribute(diagnostics, "total_ifp_ancestry_incomplete_events",
    total_ifp_ancestry_incomplete_events_);
  write_attribute(diagnostics,
    "total_ifp_ancestry_invalid_delayed_group_events",
    total_ifp_ancestry_invalid_delayed_group_events_);
  write_attribute(diagnostics, "active_batches_scored",
    static_cast<int>(batches_.size()));
  write_dataset(diagnostics, "numerator",
    std::vector<double>(result.numerators.begin(), result.numerators.end()));
  write_dataset(diagnostics, "denominator",
    std::vector<double> {result.denominator});
  H5Gclose(diagnostics);

  hid_t method_group = create_group(file_id, "method");
  write_method_group(method_group, "fclutch_cell_material", result, false);
  write_method_group(method_group, "cclutch", cclutch_result, true);
  write_clutch_ifp_method_group(method_group, clutch_ifp_result);
  write_ifp_ancestry_method_group(method_group, ifp_ancestry_result);
  H5Gclose(method_group);
  file_close(file_id);

  std::cout << "\n  beta_eff results (cell-material I*):" << std::endl;
  std::cout << "  " << std::string(136, '-') << std::endl;
  std::cout << "    group      F-CLUTCH      unc_F        C-CLUTCH      unc_C        CLUTCH-IFP    unc_CI       IFP-ancestry  unc_IFP"
            << std::endl;
  std::cout << "  " << std::string(136, '-') << std::endl;
  for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
    std::cout << "      " << std::setw(2) << (k + 1) << "      "
              << std::scientific << std::setprecision(5) << result.beta_i[k]
              << "    " << result.uncertainty[k] << "    "
              << cclutch_result.beta_i[k] << "    "
              << cclutch_result.uncertainty[k] << "    "
              << clutch_ifp_result.beta_i[k] << "    "
              << clutch_ifp_result.uncertainty[k] << "    "
              << ifp_ancestry_result.beta_i[k] << "    "
              << ifp_ancestry_result.uncertainty[k] << std::endl;
  }
  std::cout << "  " << std::string(136, '-') << std::endl;
  std::cout << "    total   " << std::scientific << std::setprecision(5)
            << result.beta_total << "    " << result.beta_total_uncertainty
            << "    " << cclutch_result.beta_total << "    "
            << cclutch_result.beta_total_uncertainty << "    "
            << clutch_ifp_result.beta_total << "    "
            << clutch_ifp_result.beta_total_uncertainty << "    "
            << ifp_ancestry_result.beta_total << "    "
            << ifp_ancestry_result.beta_total_uncertainty
            << std::endl;
  std::cout << "  " << std::string(136, '-') << std::endl;
  std::cout << "  (*) Primary method: F-CLUTCH fission-site scoring with "
               "material-resolved I*(cell,material) and same-cell missing-I "
               "fallback; C-CLUTCH folds active transfer functions with the "
               "same source-state importance; CLUTCH-IFP replaces I* by an "
               "IFP-derived source-state importance; IFP-ancestry is a direct "
               "N-generation diagnostic estimator."
            << std::endl;
}

void BetaEffectiveAccumulator::write_generation_time_to_file(
  const std::string& filename) const
{
  if (batch_active_) {
    fatal_error("BetaEffectiveAccumulator: cannot write generation time while "
                "a batch is active.");
  }

  auto result = compute_generation_time_result();
  auto cclutch_result = compute_cclutch_generation_time_result();
  if (!result.available) {
    fatal_error("F-CLUTCH generation-time denominator is zero. Verify the "
                "spatial adjoint source and active fission-site scoring.");
  }
  if (!cclutch_result.available) {
    fatal_error("C-CLUTCH generation-time denominator is zero. Verify active "
                "source birth recording and C-CLUTCH fission event scoring.");
  }

  hid_t file_id = file_open(filename, 'w');
  write_attribute(file_id, "filetype", "generation_time");
  write_attribute(file_id, "version", "2.0");
  write_attribute(file_id, "description",
    "CLUTCH-weighted neutron generation time estimates");
  write_attribute(file_id, "primary_method", "fclutch_cell_material");
  write_attribute(file_id, "units", "s");
  write_attribute(file_id, "uses_ifp", 0);
  write_attribute(file_id, "delay_treatment", "group_mean_1_over_lambda");
  write_attribute(file_id, "source_state_definition", "cell_material");
  write_attribute(file_id, "state_indexing",
    "state=cell*n_materials+material_index");
  write_attribute(file_id, "batchwise_ratio_uncertainty", 1);

  write_dataset(
    file_id, "generation_time", std::vector<double> {result.transport_lifetime});
  write_dataset(file_id, "generation_time_uncertainty",
    std::vector<double> {result.transport_lifetime_uncertainty});
  write_dataset(file_id, "emission_adjusted_generation_time",
    std::vector<double> {result.emission_adjusted_lifetime});
  write_dataset(file_id, "emission_adjusted_generation_time_uncertainty",
    std::vector<double> {result.emission_adjusted_lifetime_uncertainty});

  hid_t metadata = create_group(file_id, "metadata");
  write_attribute(metadata, "adjoint_source",
    "fission_matrix.h5/adjoint_source generated by I*=(1/k)F^T I* on the "
    "energy-integrated cell-material fission matrix");
  write_attribute(metadata, "n_materials", n_materials_);
  write_attribute(
    metadata, "n_source_states", static_cast<int64_t>(n_source_states_));
  write_dataset(metadata, "grid_shape", grid_->shape());
  write_dataset(metadata, "grid_lower_left", grid_->origin());
  write_dataset(metadata, "grid_upper_right", grid_->upper_bound());
  std::array<double, 3> pitch {grid_->pitch(), grid_->pitch(), grid_->pitch()};
  write_dataset(metadata, "grid_pitch", pitch);
  std::vector<int> material_ids;
  material_ids.reserve(n_materials_);
  for (int m = 0; m < n_materials_; ++m) {
    if (m < static_cast<int>(model::materials.size()) && model::materials[m]) {
      material_ids.push_back(model::materials[m]->id());
    } else {
      material_ids.push_back(m);
    }
  }
  write_dataset(metadata, "material_ids", material_ids);
  H5Gclose(metadata);

  hid_t method_group = create_group(file_id, "method");
  write_generation_time_method_group(
    method_group, "fclutch_cell_material", result, false);
  write_generation_time_method_group(
    method_group, "cclutch", cclutch_result, true);
  H5Gclose(method_group);
  file_close(file_id);

  std::cout << "\n  neutron generation time results (CLUTCH-weighted):"
            << std::endl;
  std::cout << "  " << std::string(96, '-') << std::endl;
  std::cout << "    method          transport_lifetime    unc_transport    "
               "emission_adjusted    unc_adjusted"
            << std::endl;
  std::cout << "  " << std::string(96, '-') << std::endl;
  std::cout << "    F-CLUTCH        " << std::scientific
            << std::setprecision(5) << result.transport_lifetime << "        "
            << result.transport_lifetime_uncertainty << "        "
            << result.emission_adjusted_lifetime << "        "
            << result.emission_adjusted_lifetime_uncertainty << std::endl;
  std::cout << "    C-CLUTCH        " << cclutch_result.transport_lifetime
            << "        " << cclutch_result.transport_lifetime_uncertainty
            << "        " << cclutch_result.emission_adjusted_lifetime
            << "        "
            << cclutch_result.emission_adjusted_lifetime_uncertainty
            << std::endl;
  std::cout << "  " << std::string(96, '-') << std::endl;
  std::cout << "  (*) transport_lifetime uses Particle::lifetime(); "
               "emission_adjusted adds group-mean delayed precursor delay "
               "1/lambda_g without sampling extra random numbers."
            << std::endl;
}

} // namespace openmc
