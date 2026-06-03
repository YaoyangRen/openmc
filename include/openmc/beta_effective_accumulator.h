#ifndef OPENMC_BETA_EFFECTIVE_ACCUMULATOR_H
#define OPENMC_BETA_EFFECTIVE_ACCUMULATOR_H

#include "hdf5.h"
#include "openmc/position.h"

#include <array>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace openmc {

class SharedMeshGrid;

class BetaEffectiveAccumulator {
public:
  static constexpr int N_DELAYED_GROUPS = 8;

  struct MethodResult {
    bool available {false};
    double denominator {0.0};
    std::array<double, N_DELAYED_GROUPS> numerators {};
    std::array<double, N_DELAYED_GROUPS> beta_i {};
    std::array<double, N_DELAYED_GROUPS> uncertainty {};
    double beta_total {0.0};
    double beta_total_uncertainty {0.0};
  };

  struct GenerationTimeResult {
    bool available {false};
    double denominator {0.0};
    double lifetime_numerator {0.0};
    double emission_adjusted_lifetime_numerator {0.0};
    double transport_lifetime {0.0};
    double transport_lifetime_uncertainty {0.0};
    double emission_adjusted_lifetime {0.0};
    double emission_adjusted_lifetime_uncertainty {0.0};
  };

  explicit BetaEffectiveAccumulator(std::shared_ptr<SharedMeshGrid> grid);

  void set_adjoint_source_spatial(
    const std::vector<double>& adjoint_source_spatial);

  void begin_batch(int batch_id);
  void end_batch(int batch_id);

  void record_source_birth(const Position& r, int64_t source_particle_id);

  void score_fission_site(const Position& r, double site_weight,
    int delayed_group, double neutron_lifetime, double delayed_group_delay);

  void score_cclutch_fission_event(const Position& r,
    int64_t source_particle_id, double total_contribution,
    const std::array<double, N_DELAYED_GROUPS>& delayed_contributions,
    double neutron_lifetime,
    const std::array<double, N_DELAYED_GROUPS>& delayed_group_delays);

  MethodResult compute_result() const;
  MethodResult compute_cclutch_result() const;
  GenerationTimeResult compute_generation_time_result() const;
  GenerationTimeResult compute_cclutch_generation_time_result() const;

  size_t n_batches() const;
  size_t scored_sites() const;
  bool has_adjoint_source() const { return source_ready_; }

  void write_to_file(const std::string& filename = "beta_eff.h5") const;
  void write_generation_time_to_file(
    const std::string& filename = "generation_time.h5") const;

private:
  struct BatchScore {
    int batch_id {-1};
    double denominator {0.0};
    std::array<double, N_DELAYED_GROUPS> numerator {};
    double lifetime_numerator {0.0};
    double emission_adjusted_lifetime_numerator {0.0};
    double cclutch_denominator {0.0};
    std::array<double, N_DELAYED_GROUPS> cclutch_numerator {};
    double cclutch_lifetime_numerator {0.0};
    double cclutch_emission_adjusted_lifetime_numerator {0.0};
    int64_t fission_sites {0};
    int64_t scored_sites {0};
    int64_t dropped_sites {0};
    int64_t invalid_delayed_group_sites {0};
    int64_t cclutch_events {0};
    int64_t cclutch_scored_events {0};
    int64_t cclutch_dropped_events {0};
    int64_t cclutch_missing_source_events {0};
  };

  int position_to_index(const Position& r) const;
  MethodResult compute_result(bool use_cclutch) const;
  GenerationTimeResult compute_generation_time_result(bool use_cclutch) const;
  void fold_current_cclutch_batch();
  void write_method_group(
    hid_t parent, const char* name, const MethodResult& result,
    bool use_cclutch) const;
  void write_generation_time_method_group(hid_t parent, const char* name,
    const GenerationTimeResult& result, bool use_cclutch) const;
  void write_batch_matrix(hid_t group, const char* name,
    const std::vector<double>& flat, hsize_t n_batches) const;

  std::shared_ptr<SharedMeshGrid> grid_;
  std::array<double, 3> upper_bound_ {};
  double inv_pitch_ {1.0};

  std::vector<double> adjoint_source_spatial_;
  bool source_ready_ {false};

  mutable std::mutex mutex_;
  bool batch_active_ {false};
  BatchScore current_batch_;
  std::vector<BatchScore> batches_;

  std::unordered_map<int64_t, int> current_source_cells_;
  std::vector<int> current_source_counts_;
  std::vector<double> current_cclutch_transfer_total_;
  std::vector<double> current_cclutch_transfer_lifetime_;
  std::vector<double> current_cclutch_transfer_emission_adjusted_lifetime_;
  std::vector<std::array<double, N_DELAYED_GROUPS>>
    current_cclutch_transfer_delayed_;

  int64_t total_fission_sites_ {0};
  int64_t total_scored_sites_ {0};
  int64_t total_dropped_sites_ {0};
  int64_t total_invalid_delayed_group_sites_ {0};
  int64_t total_cclutch_events_ {0};
  int64_t total_cclutch_scored_events_ {0};
  int64_t total_cclutch_dropped_events_ {0};
  int64_t total_cclutch_missing_source_events_ {0};
};

} // namespace openmc

#endif // OPENMC_BETA_EFFECTIVE_ACCUMULATOR_H
