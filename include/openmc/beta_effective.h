#ifndef OPENMC_BETA_EFFECTIVE_H
#define OPENMC_BETA_EFFECTIVE_H

#include "openmc/material_nuclear_data.h"

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace openmc {

struct MethodResult {
  std::array<double, N_DELAYED_GROUPS> beta_i {};
  std::array<double, N_DELAYED_GROUPS> numerators {};
  double denominator {0.0};
  double beta_total {0.0};
  bool available {false};
};

struct Position;

enum class WeightingMode {
  VOLUME_WEIGHTED,
  REACTION_RATE_WEIGHTED
};

// Computes beta_eff with fission neutron birth-energy source importance:
//
// D = sum_c dV sum_gin sum_gb phi(c,gin) I*(c,gb)
//     [P_prompt(c,gin,gb) + sum_k P_delayed(c,gin,k,gb)]
// N_k = sum_c dV sum_gin sum_gb phi(c,gin) I*(c,gb)
//     P_delayed(c,gin,k,gb)
//
// I*(c,gb) is read from fission_matrix.h5/adjoint_source_grouped.
class BetaEffective {
public:
  explicit BetaEffective(int n_sample_points = 27,
    WeightingMode weighting_mode = WeightingMode::REACTION_RATE_WEIGHTED,
    double ref_energy = 2.0e6, double ref_temperature = 293.6)
    : n_sample_points_(n_sample_points), weighting_mode_(weighting_mode),
      ref_energy_(ref_energy), ref_temperature_(ref_temperature)
  {}

  void compute_from_files(const std::string& flux_file,
    const std::string& adjoint_flux_file, const std::string& output_file,
    const std::string& fission_matrix_file = "fission_matrix.h5");

  double get_beta_i(int group) const;

  double get_beta_total() const { return beta_total_; }

  const std::array<double, N_DELAYED_GROUPS>& get_all_beta_i() const
  {
    return beta_i_;
  }

private:
  double compute_denominator_birth_spectrum(
    const std::unordered_map<int, double>& flux, double volume);

  double compute_delayed_numerator_birth_spectrum(int group,
    const std::unordered_map<int, double>& flux, double volume);

  MaterialNuclearData extract_material_nuclear_data(int material_id) const;

  void build_cell_material_map(const std::unordered_map<int, double>& flux);

  Position sample_cell_point(const Position& lower, const Position& upper,
    int point_idx, int n_points) const;

  MaterialNuclearData compute_weighted_nuclear_data(
    const std::unordered_map<int, int>& material_counts, int total_samples,
    const std::unordered_map<int, MaterialNuclearData>& material_cache) const;

  double kahan_sum(const std::vector<double>& values) const;

  std::array<double, 3> mesh_index_to_position(int mesh_idx) const;

  std::array<int, 3> get_grid_indices(int cell_idx) const;

  Position grid_index_to_position(int cell_idx) const;

  void read_flux_data(const std::string& filename,
    std::unordered_map<int, double>& flux_map, std::array<int, 3>& shape,
    double& pitch);

  void read_adjoint_flux_data(const std::string& filename,
    std::unordered_map<int, double>& adjoint_flux_map,
    std::array<int, 3>& shape, double& pitch);

  void read_fission_adjoint_source_data(const std::string& filename,
    const std::array<int, 3>& expected_shape, double expected_pitch);

  void validate_group_metadata();

  void initialize_flux_spectrum_weights();

  int group_index_from_energy(double energy_eV) const;

  size_t delayed_offset(int energy_group, int delayed_group) const;

  void write_to_file(const std::string& filename) const;

  int n_sample_points_;
  WeightingMode weighting_mode_;
  double ref_energy_;
  double ref_temperature_;

  std::array<double, N_DELAYED_GROUPS> beta_i_;
  double beta_total_;
  std::array<double, N_DELAYED_GROUPS> numerators_;
  double denominator_;

  int flux_n_groups_ {1};
  int adjoint_n_groups_ {1};
  int n_energy_groups_ {1};
  std::vector<double> flux_energy_edges_;
  std::vector<double> adjoint_energy_edges_;
  std::vector<double> energy_edges_common_;
  std::unordered_map<int, std::vector<double>> flux_group_map_;
  std::unordered_map<int, std::vector<double>> adjoint_group_map_;
  bool flux_has_group_data_ {false};
  bool adjoint_has_group_data_ {false};
  std::vector<double> flux_collapse_weights_;

  int n_heterogeneous_cells_ {0};

  MethodResult result_birth_spectrum_;

  std::array<std::vector<double>, N_DELAYED_GROUPS>
    numerator_by_group_energy_;
  std::vector<double> denominator_by_group_energy_;

  bool has_birth_adjoint_source_ {false};
  int birth_source_n_groups_ {1};
  std::vector<double> birth_source_energy_edges_;
  std::vector<double> birth_adjoint_source_grouped_;
  std::string fission_matrix_file_ {"fission_matrix.h5"};

  std::array<int, 3> grid_shape_;
  double grid_pitch_;
  std::array<double, 3> grid_pitch_vector_;
  std::array<double, 3> grid_lower_left_;
  std::array<double, 3> grid_upper_right_;

  std::unordered_map<int, MaterialNuclearData> cell_nuclear_data_;
  std::unordered_map<int, int> cell_to_material_;
  std::vector<MaterialNuclearData> unique_materials_;
};

} // namespace openmc

#endif // OPENMC_BETA_EFFECTIVE_H
