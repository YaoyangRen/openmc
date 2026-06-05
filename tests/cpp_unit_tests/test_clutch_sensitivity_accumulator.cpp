#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <unordered_map>

#include "openmc/clutch_sensitivity_accumulator.h"
#include "openmc/mesh_init.h"

using namespace openmc;

TEST_CASE("CLUTCH sensitivity accumulator scores fixed response terms")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, true);

  ClutchSensitivityAccumulator acc(
    grid, {101, 102}, {2.0, 4.0}, "hybrid");

  std::unordered_map<int64_t, double> spatial;
  spatial[0] = 10.0;
  spatial[4] = 20.0;
  acc.set_adjoint_source_spatial(spatial);

  acc.begin_batch(1);
  acc.score_contribution({0.25, 0.25, 0.25}, 1.0, {0.5, -0.25},
    {0.2, -0.1}, {0.3, -0.15}, {0.1, 0.05});
  acc.end_batch(1);

  acc.begin_batch(2);
  acc.score_contribution({1.25, 0.25, 0.25}, 2.0, {-0.1, 0.2},
    {-0.04, 0.08}, {-0.06, 0.12}, {0.05, -0.1});
  acc.end_batch(2);

  auto result = acc.compute_result();
  REQUIRE(result.available);
  REQUIRE_THAT(result.denominator,
    Catch::Matchers::WithinRel(25.0, 1.0e-12));

  REQUIRE_THAT(result.numerator[0],
    Catch::Matchers::WithinRel(2.0, 1.0e-12));
  REQUIRE_THAT(result.numerator[1],
    Catch::Matchers::WithinRel(1.0, 1.0e-12));
  REQUIRE_THAT(result.dlogk_dparameter[0],
    Catch::Matchers::WithinRel(0.08, 1.0e-12));
  REQUIRE_THAT(result.dlogk_dparameter[1],
    Catch::Matchers::WithinRel(0.04, 1.0e-12));
  REQUIRE_THAT(result.sensitivity[0],
    Catch::Matchers::WithinRel(0.16, 1.0e-12));
  REQUIRE_THAT(result.sensitivity[1],
    Catch::Matchers::WithinRel(0.16, 1.0e-12));

  REQUIRE_THAT(result.uncertainty[0],
    Catch::Matchers::WithinRel(0.416, 1.0e-12));
  REQUIRE_THAT(result.uncertainty[1],
    Catch::Matchers::WithinRel(0.384, 1.0e-12));
}

TEST_CASE("C-CLUTCH sensitivity accumulator folds transfer response terms")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, true);

  ClutchSensitivityAccumulator acc(
    grid, {101, 102}, {2.0, 4.0}, "hybrid");

  std::unordered_map<int64_t, double> spatial;
  spatial[0] = 10.0;
  spatial[4] = 20.0;
  acc.set_adjoint_source_spatial(spatial);

  acc.begin_batch(1);
  acc.record_source_birth({0.25, 0.25, 0.25}, 101);
  acc.score_cclutch_contribution(
    {0.25, 0.25, 0.25}, 101, 3.0, {0.5, -0.25});
  acc.end_batch(1);

  acc.begin_batch(2);
  acc.record_source_birth({1.25, 0.25, 0.25}, 201);
  acc.record_source_birth({1.25, 0.25, 0.25}, 202);
  acc.score_cclutch_contribution(
    {1.25, 0.25, 0.25}, 201, 4.0, {-0.1, 0.2});
  acc.score_cclutch_contribution(
    {1.25, 0.25, 0.25}, 202, 2.0, {0.3, -0.4});
  acc.end_batch(2);

  auto result = acc.compute_cclutch_result();
  REQUIRE(result.available);
  REQUIRE_THAT(result.denominator,
    Catch::Matchers::WithinRel(45.0, 1.0e-12));
  REQUIRE_THAT(result.numerator[0],
    Catch::Matchers::WithinRel(8.5, 1.0e-12));
  REQUIRE_THAT(result.numerator[1],
    Catch::Matchers::WithinRel(-3.75, 1.0e-12));
  REQUIRE_THAT(result.dlogk_dparameter[0],
    Catch::Matchers::WithinRel(8.5 / 45.0, 1.0e-12));
  REQUIRE_THAT(result.dlogk_dparameter[1],
    Catch::Matchers::WithinRel(-3.75 / 45.0, 1.0e-12));
  REQUIRE_THAT(result.sensitivity[0],
    Catch::Matchers::WithinRel(17.0 / 45.0, 1.0e-12));
  REQUIRE_THAT(result.sensitivity[1],
    Catch::Matchers::WithinRel(-15.0 / 45.0, 1.0e-12));
}
