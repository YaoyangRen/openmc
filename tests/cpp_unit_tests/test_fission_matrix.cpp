#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "openmc/fission_matrix.h"
#include "openmc/mesh_init.h"

#include <cstdio>

using namespace openmc;

TEST_CASE("Fission matrix skips early inactive batches")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, true);

  FissionMatrix matrix(grid, 5, {}, 3, 2);

  matrix.start_new_batch(1);
  matrix.record_source_birth({0.25, 0.25, 0.25}, 1, 0, 1.0);
  matrix.record_fission_site({0.25, 0.25, 0.25}, 1.0, 1, 0);

  matrix.start_new_batch(2);
  matrix.record_source_birth({0.25, 0.25, 0.25}, 2, 0, 1.0);
  matrix.record_fission_site({0.25, 0.25, 0.25}, 1.0, 2, 0);

  matrix.start_new_batch(3);
  matrix.record_source_birth({0.25, 0.25, 0.25}, 3, 1, 2.0);
  matrix.record_fission_site({0.25, 0.25, 0.25}, 2.0, 3, 1);

  matrix.start_new_batch(-1);

  REQUIRE(matrix.score_start_batch() == 3);
  REQUIRE(matrix.skipped_batches() == 2);
  REQUIRE(matrix.n_realizations() == 1);
  REQUIRE(matrix.nnz() == 1);
  REQUIRE(matrix.get_source_nonzero_cells() == 1);
  REQUIRE(matrix.get_child_nonzero_cells() == 1);
  REQUIRE(matrix.get_source_nonzero_states() == 1);
  REQUIRE(matrix.get_child_nonzero_states() == 1);
  REQUIRE(matrix.n_source_groups() == 2);
  REQUIRE(matrix.n_source_states() == static_cast<int>(grid->n_cells() * 2));

  matrix.compute_adjoint_source("uniform", 20, 1.0e-12);

  REQUIRE(matrix.is_adjoint_computed());
  REQUIRE(matrix.is_adjoint_converged());
  REQUIRE(matrix.get_adjoint_nonzero_cells() == 1);
  REQUIRE_THAT(matrix.get_adjoint_final_residual(),
    Catch::Matchers::WithinAbs(0.0, 1.0e-12));

  matrix.finalize("test_fission_matrix.h5");
  std::remove("test_fission_matrix.h5");
}

TEST_CASE("Fission matrix separates material source states")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, true);

  FissionMatrix matrix(grid, 3, {}, 1, 2);

  matrix.start_new_batch(1);
  matrix.record_source_birth({0.25, 0.25, 0.25}, 1, 0, 1.0);
  matrix.record_fission_site({0.25, 0.25, 0.25}, 1.0, 1, 1);
  matrix.record_source_birth({0.25, 0.25, 0.25}, 2, 1, 1.0);
  matrix.record_fission_site({0.25, 0.25, 0.25}, 1.0, 2, 0);
  matrix.start_new_batch(-1);

  REQUIRE(matrix.nnz() == 2);
  REQUIRE(matrix.get_source_nonzero_cells() == 1);
  REQUIRE(matrix.get_child_nonzero_cells() == 1);
  REQUIRE(matrix.get_source_nonzero_states() == 2);
  REQUIRE(matrix.get_child_nonzero_states() == 2);
  REQUIRE(matrix.n_source_states() == grid->n_cells() * 2);
  REQUIRE(matrix.get_adjoint_source().empty());
}
