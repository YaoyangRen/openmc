#ifndef OPENMC_CLUTCH_SENSITIVITY_ACCUMULATOR_H
#define OPENMC_CLUTCH_SENSITIVITY_ACCUMULATOR_H

#include "hdf5.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/vector.h"

#include <array>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace openmc {

class SharedMeshGrid;

class ClutchSensitivityAccumulator {
public:
  struct BatchScore {
    int batch_id {-1};
    double denominator {0.0};
    vector<double> numerator;
    double cclutch_denominator {0.0};
    vector<double> cclutch_numerator;
    int64_t fission_sites {0};
    int64_t scored_sites {0};
    int64_t dropped_sites {0};
    int64_t cclutch_events {0};
    int64_t cclutch_scored_events {0};
    int64_t cclutch_dropped_events {0};
    int64_t cclutch_missing_source_events {0};
  };

  struct MethodResult {
    bool available {false};
    double denominator {0.0};
    vector<double> numerator;
    vector<double> dlogk_dparameter;
    vector<double> sensitivity;
    vector<double> uncertainty;
  };

  ClutchSensitivityAccumulator(std::shared_ptr<SharedMeshGrid> grid,
    vector<int> derivative_indices, std::string method);
  ClutchSensitivityAccumulator(std::shared_ptr<SharedMeshGrid> grid,
    vector<int> parameter_ids, vector<double> parameter_values,
    std::string method);

  void set_adjoint_source_spatial(
    const vector<double>& adjoint_source_spatial);

  void begin_batch(int batch_id);
  void end_batch(int batch_id);

  void record_source_birth(const Position& r, int64_t source_particle_id);

  void score_fission_site(
    const Particle& p, const Position& r, double site_weight);

  void score_cclutch_fission_event(const Particle& p, const Position& r,
    double total_contribution);

  void score_cclutch_contribution(const Position& r,
    int64_t source_particle_id, double total_contribution,
    const vector<double>& response_terms);

  void score_contribution(const Position& r, double site_weight,
    const vector<double>& path_terms, const vector<double>& track_terms,
    const vector<double>& collision_terms, const vector<double>& direct_terms);

  MethodResult compute_result(const vector<BatchScore>& batches) const;
  MethodResult compute_cclutch_result(const vector<BatchScore>& batches) const;
  MethodResult compute_result() const;
  MethodResult compute_cclutch_result() const;
  void write_to_file(const std::string& filename) const;

  size_t n_parameters() const { return derivative_indices_.size(); }
  bool has_adjoint_source() const { return source_ready_; }

private:
  int position_to_index(const Position& r) const;
  vector<BatchScore> collect_global_batches(
    int64_t& total_fission_sites, int64_t& total_scored_sites,
    int64_t& total_dropped_sites, int64_t& total_cclutch_events,
    int64_t& total_cclutch_scored_events,
    int64_t& total_cclutch_dropped_events,
    int64_t& total_cclutch_missing_source_events) const;
  void write_batch_matrix(hid_t group, const char* name,
    const vector<double>& flat, hsize_t n_batches) const;
  void write_method_group(hid_t parent, const char* name,
    const MethodResult& result, const vector<BatchScore>& batches,
    bool use_cclutch) const;
  void initialize_parameter_metadata();
  void fold_current_cclutch_batch();

  std::shared_ptr<SharedMeshGrid> grid_;
  std::array<double, 3> upper_bound_ {};
  double inv_pitch_ {1.0};

  vector<double> adjoint_source_spatial_;
  bool source_ready_ {false};

  vector<int> derivative_indices_;
  vector<int> derivative_ids_;
  vector<std::string> derivative_variables_;
  vector<int> derivative_material_ids_;
  vector<std::string> derivative_nuclides_;
  vector<double> parameter_values_;
  std::string method_;

  mutable std::mutex mutex_;
  bool batch_active_ {false};
  BatchScore current_batch_;
  vector<BatchScore> batches_;

  std::unordered_map<int64_t, int> current_source_cells_;
  vector<int> current_source_counts_;
  vector<double> current_cclutch_transfer_total_;
  vector<double> current_cclutch_transfer_numerator_;

  int64_t total_fission_sites_ {0};
  int64_t total_scored_sites_ {0};
  int64_t total_dropped_sites_ {0};
  int64_t total_cclutch_events_ {0};
  int64_t total_cclutch_scored_events_ {0};
  int64_t total_cclutch_dropped_events_ {0};
  int64_t total_cclutch_missing_source_events_ {0};
};

} // namespace openmc

#endif // OPENMC_CLUTCH_SENSITIVITY_ACCUMULATOR_H
