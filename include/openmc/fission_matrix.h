#ifndef OPENMC_FISSION_MATRIX_H
#define OPENMC_FISSION_MATRIX_H

#include "openmc/mesh_init.h"
#include "openmc/position.h"
#include "openmc/vector.h"

#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>

namespace openmc {

class FissionMatrix {
public:
  explicit FissionMatrix(std::shared_ptr<SharedMeshGrid> grid, int max_batches,
    std::vector<double> energy_edges = {});

  void record_source_birth(const Position& r, int64_t source_particle_id,
    double energy = -1.0, int mg_group = -1, double source_weight = 1.0);

  void record_fission_site(const Position& r, double source_weight,
    int64_t source_particle_id, double energy = -1.0, int mg_group = -1);

  void start_new_batch(int batch_id);
  void finalize(const std::string& filename = "fission_matrix.h5");

  void compute_adjoint_source(const std::string& initial_guess = "uniform",
    int max_iterations = 1, double tolerance = 1.0e-6);

  bool is_adjoint_computed() const { return adjoint_computed_; }
  int get_adjoint_nonzero_cells() const;

  const vector<double>& get_adjoint_source() const { return adjoint_source_; }
  const vector<double>& get_spatial_adjoint_source() const
  {
    return adjoint_source_;
  }

  int n_source_groups() const { return 1; }

  const std::array<int, 3>& shape() const { return grid_->shape(); }
  const std::array<double, 3>& origin() const { return grid_->origin(); }
  double pitch() const { return grid_->pitch(); }
  size_t n_cells() const { return grid_->n_cells(); }

private:
  int position_to_index(const Position& r) const;

  std::shared_ptr<SharedMeshGrid> grid_;
  std::array<double, 3> upper_bound_ {};
  double inv_pitch_ {1.0};

  int current_batch_id_ {-1};
  int max_batches_ {0};
  int n_realizations_ {0};

  // key = parent_cell * n_cells + child_cell
  std::unordered_map<size_t, double> current_batch_sparse_;
  std::unordered_map<size_t, double> fission_matrix_sparse_;
  std::unordered_map<int64_t, int> source_birth_cells_;

  vector<double> source_counts_;
  vector<double> current_batch_source_counts_;
  vector<double> adjoint_source_;

  std::atomic<uint64_t> total_fissions_ {0};
  std::atomic<uint64_t> total_sources_ {0};

  bool adjoint_computed_ {false};
  int adjoint_iterations_ {0};
  double keff_reference_ {1.0};

  mutable std::mutex data_mutex_;
};

} // namespace openmc

#endif // OPENMC_FISSION_MATRIX_H
