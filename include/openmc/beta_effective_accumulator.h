#ifndef OPENMC_BETA_EFFECTIVE_ACCUMULATOR_H
#define OPENMC_BETA_EFFECTIVE_ACCUMULATOR_H

#include "hdf5.h"
#include "openmc/position.h"

#include <array>
#include <memory>
#include <mutex>
#include <string>
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

  explicit BetaEffectiveAccumulator(std::shared_ptr<SharedMeshGrid> grid);

  void set_adjoint_source_spatial(
    const std::vector<double>& adjoint_source_spatial);

  void begin_batch(int batch_id);
  void end_batch(int batch_id);

  void score_fission_site(
    const Position& r, double site_weight, int delayed_group);

  MethodResult compute_result() const;

  size_t n_batches() const;
  size_t scored_sites() const;
  bool has_adjoint_source() const { return source_ready_; }

  void write_to_file(const std::string& filename = "beta_eff.h5") const;

private:
  struct BatchScore {
    int batch_id {-1};
    double denominator {0.0};
    std::array<double, N_DELAYED_GROUPS> numerator {};
    int64_t fission_sites {0};
    int64_t scored_sites {0};
    int64_t dropped_sites {0};
    int64_t invalid_delayed_group_sites {0};
  };

  int position_to_index(const Position& r) const;
  void write_method_group(hid_t parent, const MethodResult& result) const;
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

  int64_t total_fission_sites_ {0};
  int64_t total_scored_sites_ {0};
  int64_t total_dropped_sites_ {0};
  int64_t total_invalid_delayed_group_sites_ {0};
};

} // namespace openmc

#endif // OPENMC_BETA_EFFECTIVE_ACCUMULATOR_H
