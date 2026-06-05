#ifndef OPENMC_FISSION_MATRIX_H
#define OPENMC_FISSION_MATRIX_H

#include "openmc/mesh_init.h"
#include "openmc/position.h"
#include "openmc/vector.h"

#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace openmc {

class FissionMatrix {
public:
  explicit FissionMatrix(std::shared_ptr<SharedMeshGrid> grid, int max_batches,
    std::vector<double> energy_edges = {}, int score_start_batch = 1,
    int n_materials = -1);

  void record_source_birth(const Position& r, int64_t source_particle_id,
    int material_index, double source_weight = 1.0);

  void record_fission_site(const Position& r, double source_weight,
    int64_t source_particle_id, int material_index);

  void start_new_batch(int batch_id);
  void finalize(const std::string& filename = "fission_matrix.h5");

  void compute_adjoint_source(const std::string& initial_guess = "uniform",
    int max_iterations = 1, double tolerance = 1.0e-6);

  bool is_adjoint_computed() const { return adjoint_computed_; }
  bool is_adjoint_converged() const { return adjoint_converged_; }
  int get_adjoint_nonzero_cells() const;
  int get_adjoint_nonzero_states() const;
  int get_source_nonzero_cells() const;
  int get_child_nonzero_cells() const;
  int get_source_nonzero_states() const;
  int get_child_nonzero_states() const;
  int score_start_batch() const { return score_start_batch_; }
  int n_realizations() const { return n_realizations_; }
  int skipped_batches() const { return skipped_batches_; }
  size_t nnz() const { return fission_matrix_sparse_.size(); }
  double get_adjoint_final_residual() const
  {
    return adjoint_final_residual_;
  }

  const std::unordered_map<int64_t, double>& get_adjoint_source() const
  {
    return adjoint_source_;
  }
  const std::unordered_map<int64_t, double>& get_spatial_adjoint_source() const
  {
    return adjoint_source_;
  }
  std::unordered_map<int64_t, double> get_cell_adjoint_source() const;

  int n_source_groups() const { return n_materials_; }
  size_t n_source_states() const { return n_source_states_; }

  const std::array<int, 3>& shape() const { return grid_->shape(); }
  const std::array<double, 3>& origin() const { return grid_->origin(); }
  double pitch() const { return grid_->pitch(); }
  size_t n_cells() const { return grid_->n_cells(); }

private:
  struct StatePairKey {
    int64_t parent_state;
    int64_t child_state;

    bool operator==(const StatePairKey& other) const
    {
      return parent_state == other.parent_state &&
             child_state == other.child_state;
    }
  };

  struct StatePairKeyHash {
    size_t operator()(const StatePairKey& key) const
    {
      size_t seed = std::hash<int64_t> {}(key.parent_state);
      const size_t child_hash = std::hash<int64_t> {}(key.child_state);
      seed ^= child_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };

  int position_to_index(const Position& r) const;
  int64_t source_state_index(const Position& r, int material_index) const;
  bool should_score_batch(int batch_id) const;
  bool is_scoring_current_batch() const;

  std::shared_ptr<SharedMeshGrid> grid_;
  std::array<double, 3> upper_bound_ {};
  double inv_pitch_ {1.0};

  int current_batch_id_ {-1};
  int max_batches_ {0};
  int n_realizations_ {0};
  int skipped_batches_ {0};
  int score_start_batch_ {1};
  int n_materials_ {1};
  size_t n_source_states_ {0};

  std::unordered_map<StatePairKey, double, StatePairKeyHash>
    current_batch_sparse_;
  std::unordered_map<StatePairKey, double, StatePairKeyHash>
    fission_matrix_sparse_;
  std::unordered_map<int64_t, int64_t> source_birth_states_;

  std::unordered_map<int64_t, double> source_counts_;
  std::unordered_map<int64_t, double> current_batch_source_counts_;
  std::unordered_map<int64_t, double> adjoint_source_;

  std::atomic<uint64_t> total_fissions_ {0};
  std::atomic<uint64_t> total_sources_ {0};
  std::atomic<uint64_t> total_invalid_source_states_ {0};
  std::atomic<uint64_t> total_invalid_fission_states_ {0};

  bool adjoint_computed_ {false};
  bool adjoint_converged_ {false};
  int adjoint_iterations_ {0};
  double adjoint_final_residual_ {0.0};
  double keff_reference_ {1.0};

  mutable std::mutex data_mutex_;
};

} // namespace openmc

#endif // OPENMC_FISSION_MATRIX_H
