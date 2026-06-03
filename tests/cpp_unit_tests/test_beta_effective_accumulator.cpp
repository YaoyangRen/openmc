#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>

#include "openmc/beta_effective_accumulator.h"
#include "openmc/mesh_init.h"

using namespace openmc;

TEST_CASE("F-CLUTCH beta_eff accumulator scores spatial I-star")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, true);

  BetaEffectiveAccumulator acc(grid);

  std::vector<double> spatial(grid->n_cells(), 0.0);
  spatial[0] = 10.0;
  spatial[4] = 20.0;
  acc.set_adjoint_source_spatial(spatial);

  acc.begin_batch(1);
  acc.score_fission_site({0.25, 0.25, 0.25}, 1.0, 0);
  acc.score_fission_site({0.25, 0.25, 0.25}, 2.0, 1);
  acc.end_batch(1);

  acc.begin_batch(2);
  acc.score_fission_site({1.25, 0.25, 0.25}, 4.0, 2);
  acc.end_batch(2);

  auto result = acc.compute_result();
  REQUIRE(result.available);
  REQUIRE_THAT(result.denominator,
    Catch::Matchers::WithinRel(55.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[0],
    Catch::Matchers::WithinRel(10.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[1],
    Catch::Matchers::WithinRel(40.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[0],
    Catch::Matchers::WithinRel(10.0 / 55.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[1],
    Catch::Matchers::WithinRel(40.0 / 55.0, 1.0e-12));
  REQUIRE_THAT(result.beta_total,
    Catch::Matchers::WithinRel(50.0 / 55.0, 1.0e-12));
  REQUIRE(result.uncertainty[0] > 0.0);
  REQUIRE(result.uncertainty[1] > 0.0);
}

TEST_CASE("C-CLUTCH beta_eff accumulator folds source-normalized transfer")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, true);

  BetaEffectiveAccumulator acc(grid);

  std::vector<double> spatial(grid->n_cells(), 0.0);
  spatial[0] = 10.0;
  spatial[4] = 20.0;
  acc.set_adjoint_source_spatial(spatial);

  std::array<double, BetaEffectiveAccumulator::N_DELAYED_GROUPS> delayed {};

  acc.begin_batch(1);
  acc.record_source_birth({0.25, 0.25, 0.25}, 101);
  delayed.fill(0.0);
  delayed[0] = 0.3;
  delayed[1] = 0.6;
  acc.score_cclutch_fission_event({0.25, 0.25, 0.25}, 101, 3.0, delayed);
  acc.end_batch(1);

  acc.begin_batch(2);
  acc.record_source_birth({1.25, 0.25, 0.25}, 201);
  acc.record_source_birth({1.25, 0.25, 0.25}, 202);
  delayed.fill(0.0);
  delayed[1] = 1.0;
  acc.score_cclutch_fission_event({1.25, 0.25, 0.25}, 201, 4.0, delayed);
  delayed.fill(0.0);
  delayed[0] = 0.5;
  acc.score_cclutch_fission_event({1.25, 0.25, 0.25}, 202, 2.0, delayed);
  acc.end_batch(2);

  auto result = acc.compute_cclutch_result();
  REQUIRE(result.available);
  REQUIRE_THAT(result.denominator,
    Catch::Matchers::WithinRel(45.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[0],
    Catch::Matchers::WithinRel(4.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[1],
    Catch::Matchers::WithinRel(8.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[0],
    Catch::Matchers::WithinRel(4.0 / 45.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[1],
    Catch::Matchers::WithinRel(8.0 / 45.0, 1.0e-12));
  REQUIRE_THAT(result.beta_total,
    Catch::Matchers::WithinRel(12.0 / 45.0, 1.0e-12));
}
