#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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
