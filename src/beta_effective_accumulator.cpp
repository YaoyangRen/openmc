#include "openmc/beta_effective_accumulator.h"

#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/mesh_init.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace openmc {

BetaEffectiveAccumulator::BetaEffectiveAccumulator(
  std::shared_ptr<SharedMeshGrid> grid)
  : grid_ {std::move(grid)}
{
  if (!grid_) {
    throw std::runtime_error(
      "BetaEffectiveAccumulator requires a valid SharedMeshGrid.");
  }

  upper_bound_ = grid_->upper_bound();
  inv_pitch_ = grid_->inv_pitch();
  const size_t n_cells = grid_->n_cells();
  current_source_counts_.assign(n_cells, 0);
  current_cclutch_transfer_total_.assign(n_cells, 0.0);
  std::array<double, N_DELAYED_GROUPS> zero_delayed {};
  current_cclutch_transfer_delayed_.assign(n_cells, zero_delayed);
}

void BetaEffectiveAccumulator::set_adjoint_source_spatial(
  const std::vector<double>& adjoint_source_spatial)
{
  const size_t n_cells = grid_->n_cells();
  if (adjoint_source_spatial.size() != n_cells) {
    fatal_error("BetaEffectiveAccumulator: spatial adjoint source size "
                "does not match n_cells.");
  }

  bool has_positive = false;
  for (double v : adjoint_source_spatial) {
    if (!std::isfinite(v)) {
      fatal_error("BetaEffectiveAccumulator: spatial adjoint source contains "
                  "non-finite values.");
    }
    has_positive = has_positive || v > 0.0;
  }
  if (!has_positive) {
    fatal_error(
      "BetaEffectiveAccumulator: spatial adjoint source is all zero.");
  }

  std::lock_guard<std::mutex> lock(mutex_);
  adjoint_source_spatial_ = adjoint_source_spatial;
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
  current_source_cells_.clear();
  std::fill(current_source_counts_.begin(), current_source_counts_.end(), 0);
  std::fill(current_cclutch_transfer_total_.begin(),
    current_cclutch_transfer_total_.end(), 0.0);
  for (auto& delayed : current_cclutch_transfer_delayed_) {
    delayed.fill(0.0);
  }
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
  batches_.push_back(current_batch_);
  batch_active_ = false;
  current_batch_ = BatchScore {};
  current_source_cells_.clear();
}

void BetaEffectiveAccumulator::record_source_birth(
  const Position& r, int64_t source_particle_id)
{
  if (source_particle_id < 0) {
    return;
  }

  std::lock_guard<std::mutex> lock(mutex_);
  if (!batch_active_) {
    return;
  }

  const int cell = position_to_index(r);
  if (cell < 0) {
    return;
  }

  auto [it, inserted] = current_source_cells_.emplace(source_particle_id, cell);
  if (inserted) {
    ++current_source_counts_[cell];
  } else if (it->second != cell) {
    --current_source_counts_[it->second];
    it->second = cell;
    ++current_source_counts_[cell];
  }
}

void BetaEffectiveAccumulator::score_fission_site(
  const Position& r, double site_weight, int delayed_group)
{
  if (!std::isfinite(site_weight) || site_weight <= 0.0) {
    return;
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

  const int cell = position_to_index(r);
  if (cell < 0) {
    ++current_batch_.dropped_sites;
    ++total_dropped_sites_;
    return;
  }

  const double importance = adjoint_source_spatial_[cell];
  if (importance <= 0.0) {
    ++current_batch_.dropped_sites;
    ++total_dropped_sites_;
    return;
  }

  current_batch_.denominator += site_weight * importance;

  if (delayed_group < 0 || delayed_group > N_DELAYED_GROUPS) {
    ++current_batch_.invalid_delayed_group_sites;
    ++total_invalid_delayed_group_sites_;
  } else if (delayed_group > 0) {
    current_batch_.numerator[delayed_group - 1] += site_weight * importance;
  }

  ++current_batch_.scored_sites;
  ++total_scored_sites_;
}

void BetaEffectiveAccumulator::score_cclutch_fission_event(const Position& r,
  int64_t source_particle_id, double total_contribution,
  const std::array<double, N_DELAYED_GROUPS>& delayed_contributions)
{
  if (!std::isfinite(total_contribution) || total_contribution <= 0.0) {
    return;
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

  auto it = current_source_cells_.find(source_particle_id);
  if (it == current_source_cells_.end()) {
    ++current_batch_.cclutch_missing_source_events;
    ++total_cclutch_missing_source_events_;
    return;
  }

  const int source_cell = it->second;
  if (source_cell < 0 ||
      source_cell >= static_cast<int>(current_cclutch_transfer_total_.size())) {
    ++current_batch_.cclutch_dropped_events;
    ++total_cclutch_dropped_events_;
    return;
  }

  current_cclutch_transfer_total_[source_cell] += total_contribution;
  for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
    const double delayed = delayed_contributions[k];
    if (std::isfinite(delayed) && delayed > 0.0) {
      current_cclutch_transfer_delayed_[source_cell][k] += delayed;
    }
  }

  ++current_batch_.cclutch_scored_events;
  ++total_cclutch_scored_events_;
}

void BetaEffectiveAccumulator::fold_current_cclutch_batch()
{
  for (int source_cell = 0;
       source_cell < static_cast<int>(current_cclutch_transfer_total_.size());
       ++source_cell) {
    const double total = current_cclutch_transfer_total_[source_cell];
    bool has_delayed = false;
    for (double delayed : current_cclutch_transfer_delayed_[source_cell]) {
      has_delayed = has_delayed || delayed > 0.0;
    }
    if (total <= 0.0 && !has_delayed) {
      continue;
    }

    const int source_count = current_source_counts_[source_cell];
    if (source_count <= 0) {
      fatal_error("BetaEffectiveAccumulator: C-CLUTCH transfer response exists "
                  "for a source cell with zero source count.");
    }

    const double importance = adjoint_source_spatial_[source_cell];
    if (importance <= 0.0) {
      continue;
    }

    const double source_weight =
      importance / static_cast<double>(source_count);
    current_batch_.cclutch_denominator += source_weight * total;
    for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
      current_batch_.cclutch_numerator[k] +=
        source_weight * current_cclutch_transfer_delayed_[source_cell][k];
    }
  }
}

BetaEffectiveAccumulator::MethodResult
BetaEffectiveAccumulator::compute_result() const
{
  return compute_result(false);
}

BetaEffectiveAccumulator::MethodResult
BetaEffectiveAccumulator::compute_cclutch_result() const
{
  return compute_result(true);
}

BetaEffectiveAccumulator::MethodResult
BetaEffectiveAccumulator::compute_result(bool use_cclutch) const
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
    sum_denominator +=
      use_cclutch ? batch.cclutch_denominator : batch.denominator;
    for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
      sum_numerator[k] +=
        use_cclutch ? batch.cclutch_numerator[k] : batch.numerator[k];
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
        const double numerator =
          use_cclutch ? batch.cclutch_numerator[k] : batch.numerator[k];
        const double denominator =
          use_cclutch ? batch.cclutch_denominator : batch.denominator;
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
        use_cclutch ? batch.cclutch_numerator : batch.numerator;
      const double numerator_sum =
        std::accumulate(numerator.begin(), numerator.end(), 0.0);
      const double denominator =
        use_cclutch ? batch.cclutch_denominator : batch.denominator;
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
  write_attribute(group, "source_state_definition", "cell");
  if (use_cclutch) {
    write_attribute(group, "formula",
      "D=mean_b sum_source I*(source_cell) T_total(source_cell); "
      "N_k=mean_b sum_source I*(source_cell) T_delayed_k(source_cell)");
  } else {
    write_attribute(group, "formula",
      "D=mean_b sum_sites w_site I*(cell); "
      "N_k=mean_b sum_delayed_k w_site I*(cell)");
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

void BetaEffectiveAccumulator::write_to_file(const std::string& filename) const
{
  if (batch_active_) {
    fatal_error(
      "BetaEffectiveAccumulator: cannot write beta_eff while a batch is active.");
  }

  auto result = compute_result();
  auto cclutch_result = compute_cclutch_result();
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

  hid_t file_id = file_open(filename, 'w');
  write_attribute(file_id, "filetype", "beta_effective");
  write_attribute(file_id, "version", "3.0");
  write_attribute(
    file_id, "description", "Effective delayed neutron fraction beta_eff");

  write_dataset(file_id, "beta_i",
    std::vector<double>(result.beta_i.begin(), result.beta_i.end()));
  write_dataset(file_id, "beta_total", std::vector<double> {result.beta_total});
  write_dataset(file_id, "uncertainty",
    std::vector<double>(result.uncertainty.begin(), result.uncertainty.end()));

  hid_t metadata = create_group(file_id, "metadata");
  write_attribute(metadata, "primary_method", "fclutch_spatial");
  write_attribute(metadata, "theory_reference", "Qiu2016_F_CLUTCH_Eq31_Eq43");
  write_attribute(metadata, "uses_fission_event_sites", 1);
  write_attribute(metadata, "has_cclutch_method", 1);
  write_attribute(metadata, "uses_birth_energy_importance", 0);
  write_attribute(metadata, "source_state_definition", "cell");
  write_attribute(metadata, "n_delayed_groups", N_DELAYED_GROUPS);
  write_attribute(metadata, "batchwise_ratio_uncertainty", 1);
  write_attribute(metadata, "adjoint_source",
    "fission_matrix.h5/adjoint_source generated by I*=(1/k)F^T I* on the "
    "energy-integrated cell-to-cell fission matrix");
  write_dataset(metadata, "grid_shape", grid_->shape());
  write_dataset(metadata, "grid_lower_left", grid_->origin());
  write_dataset(metadata, "grid_upper_right", grid_->upper_bound());
  std::array<double, 3> pitch {grid_->pitch(), grid_->pitch(), grid_->pitch()};
  write_dataset(metadata, "grid_pitch", pitch);
  H5Gclose(metadata);

  hid_t diagnostics = create_group(file_id, "diagnostics");
  write_attribute(diagnostics, "total_fission_sites", total_fission_sites_);
  write_attribute(diagnostics, "total_scored_sites", total_scored_sites_);
  write_attribute(diagnostics, "total_dropped_sites", total_dropped_sites_);
  write_attribute(diagnostics, "total_invalid_delayed_group_sites",
    total_invalid_delayed_group_sites_);
  write_attribute(diagnostics, "total_cclutch_events", total_cclutch_events_);
  write_attribute(diagnostics, "total_cclutch_scored_events",
    total_cclutch_scored_events_);
  write_attribute(diagnostics, "total_cclutch_dropped_events",
    total_cclutch_dropped_events_);
  write_attribute(diagnostics, "total_cclutch_missing_source_events",
    total_cclutch_missing_source_events_);
  write_attribute(diagnostics, "active_batches_scored",
    static_cast<int>(batches_.size()));
  write_dataset(diagnostics, "numerator",
    std::vector<double>(result.numerators.begin(), result.numerators.end()));
  write_dataset(diagnostics, "denominator",
    std::vector<double> {result.denominator});
  H5Gclose(diagnostics);

  hid_t method_group = create_group(file_id, "method");
  write_method_group(method_group, "fclutch_spatial", result, false);
  write_method_group(method_group, "cclutch", cclutch_result, true);
  H5Gclose(method_group);
  file_close(file_id);

  const double mcnp_beta[6] = {
    0.00016, 0.00104, 0.00097, 0.00253, 0.00107, 0.00042};
  const double mcnp_total = 0.00621;

  std::cout << "\n  beta_eff results (spatial I*):" << std::endl;
  std::cout << "  " << std::string(112, '-') << std::endl;
  std::cout << "    group      F-CLUTCH      unc_F        C-CLUTCH      unc_C        MCNP ref     F bias(%)   C bias(%)"
            << std::endl;
  std::cout << "  " << std::string(112, '-') << std::endl;
  for (int k = 0; k < N_DELAYED_GROUPS; ++k) {
    std::cout << "      " << std::setw(2) << (k + 1) << "      "
              << std::scientific << std::setprecision(5) << result.beta_i[k]
              << "    " << result.uncertainty[k] << "    "
              << cclutch_result.beta_i[k] << "    "
              << cclutch_result.uncertainty[k];
    if (k < 6) {
      const double f_bias =
        (result.beta_i[k] - mcnp_beta[k]) / mcnp_beta[k] * 100.0;
      const double c_bias =
        (cclutch_result.beta_i[k] - mcnp_beta[k]) / mcnp_beta[k] * 100.0;
      std::cout << "    " << mcnp_beta[k] << "    " << std::fixed
                << std::setprecision(2) << std::setw(8) << f_bias << "    "
                << std::setw(8) << c_bias;
    }
    std::cout << std::endl;
  }
  const double f_total_bias =
    (result.beta_total - mcnp_total) / mcnp_total * 100.0;
  const double c_total_bias =
    (cclutch_result.beta_total - mcnp_total) / mcnp_total * 100.0;
  std::cout << "  " << std::string(112, '-') << std::endl;
  std::cout << "    total   " << std::scientific << std::setprecision(5)
            << result.beta_total << "    " << result.beta_total_uncertainty
            << "    " << cclutch_result.beta_total << "    "
            << cclutch_result.beta_total_uncertainty << "    " << mcnp_total
            << "    " << std::fixed << std::setprecision(2) << std::setw(8)
            << f_total_bias << "    " << std::setw(8) << c_total_bias
            << std::endl;
  std::cout << "  " << std::string(112, '-') << std::endl;
  std::cout << "  (*) Primary method: F-CLUTCH fission-site scoring with "
               "spatial I*(cell); C-CLUTCH folds active transfer functions "
               "with the same spatial source importance."
            << std::endl;
}

} // namespace openmc
