#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <unordered_map>

#include "openmc/bank.h"
#include "openmc/beta_effective_accumulator.h"
#include "openmc/settings.h"
#include "openmc/mesh_init.h"

using namespace openmc;

TEST_CASE("F-CLUTCH beta_eff accumulator scores spatial I-star")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, true);

  BetaEffectiveAccumulator acc(grid, 2);

  auto state = [](int cell, int material) { return cell * 2 + material; };
  std::unordered_map<int64_t, double> spatial;
  spatial[state(0, 0)] = 10.0;
  spatial[state(0, 1)] = 30.0;
  spatial[state(4, 1)] = 20.0;
  acc.set_adjoint_source_spatial(spatial);

  acc.begin_batch(1);
  acc.score_fission_site({0.25, 0.25, 0.25}, 1.0, 0, 0, 1.0, 0.0);
  acc.score_fission_site({0.25, 0.25, 0.25}, 2.0, 1, 1, 2.0, 5.0);
  acc.end_batch(1);

  acc.begin_batch(2);
  acc.score_fission_site({1.25, 0.25, 0.25}, 4.0, 1, 2, 3.0, 7.0);
  acc.end_batch(2);

  auto result = acc.compute_result();
  REQUIRE(result.available);
  REQUIRE_THAT(result.denominator,
    Catch::Matchers::WithinRel(75.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[0],
    Catch::Matchers::WithinRel(30.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[1],
    Catch::Matchers::WithinRel(40.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[0],
    Catch::Matchers::WithinRel(30.0 / 75.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[1],
    Catch::Matchers::WithinRel(40.0 / 75.0, 1.0e-12));
  REQUIRE_THAT(result.beta_total,
    Catch::Matchers::WithinRel(70.0 / 75.0, 1.0e-12));
  REQUIRE(result.uncertainty[0] > 0.0);
  REQUIRE(result.uncertainty[1] > 0.0);

  auto time = acc.compute_generation_time_result();
  REQUIRE(time.available);
  REQUIRE_THAT(time.denominator, Catch::Matchers::WithinRel(75.0, 1.0e-12));
  REQUIRE_THAT(
    time.lifetime_numerator, Catch::Matchers::WithinRel(185.0, 1.0e-12));
  REQUIRE_THAT(time.emission_adjusted_lifetime_numerator,
    Catch::Matchers::WithinRel(615.0, 1.0e-12));
  REQUIRE_THAT(time.transport_lifetime,
    Catch::Matchers::WithinRel(185.0 / 75.0, 1.0e-12));
  REQUIRE_THAT(time.emission_adjusted_lifetime,
    Catch::Matchers::WithinRel(615.0 / 75.0, 1.0e-12));
  REQUIRE(time.transport_lifetime_uncertainty > 0.0);
  REQUIRE(time.emission_adjusted_lifetime_uncertainty > 0.0);
}

TEST_CASE("C-CLUTCH beta_eff accumulator folds source-normalized transfer")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, true);

  BetaEffectiveAccumulator acc(grid, 2);

  auto state = [](int cell, int material) { return cell * 2 + material; };
  std::unordered_map<int64_t, double> spatial;
  spatial[state(0, 0)] = 10.0;
  spatial[state(4, 1)] = 20.0;
  acc.set_adjoint_source_spatial(spatial);

  std::array<double, BetaEffectiveAccumulator::N_DELAYED_GROUPS> delayed {};
  std::array<double, BetaEffectiveAccumulator::N_DELAYED_GROUPS> delays {};
  delays[0] = 5.0;
  delays[1] = 7.0;

  acc.begin_batch(1);
  acc.record_source_birth({0.25, 0.25, 0.25}, 101, 0);
  delayed.fill(0.0);
  delayed[0] = 0.3;
  delayed[1] = 0.6;
  acc.score_cclutch_fission_event(
    {0.25, 0.25, 0.25}, 101, 3.0, delayed, 2.0, delays);
  acc.end_batch(1);

  acc.begin_batch(2);
  acc.record_source_birth({1.25, 0.25, 0.25}, 201, 1);
  acc.record_source_birth({1.25, 0.25, 0.25}, 202, 1);
  delayed.fill(0.0);
  delayed[1] = 1.0;
  acc.score_cclutch_fission_event(
    {1.25, 0.25, 0.25}, 201, 4.0, delayed, 3.0, delays);
  delayed.fill(0.0);
  delayed[0] = 0.5;
  acc.score_cclutch_fission_event(
    {1.25, 0.25, 0.25}, 202, 2.0, delayed, 5.0, delays);
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

  auto time = acc.compute_cclutch_generation_time_result();
  REQUIRE(time.available);
  REQUIRE_THAT(time.denominator, Catch::Matchers::WithinRel(45.0, 1.0e-12));
  REQUIRE_THAT(
    time.lifetime_numerator, Catch::Matchers::WithinRel(140.0, 1.0e-12));
  REQUIRE_THAT(time.emission_adjusted_lifetime_numerator,
    Catch::Matchers::WithinRel(216.0, 1.0e-12));
  REQUIRE_THAT(time.transport_lifetime,
    Catch::Matchers::WithinRel(140.0 / 45.0, 1.0e-12));
  REQUIRE_THAT(time.emission_adjusted_lifetime,
    Catch::Matchers::WithinRel(216.0 / 45.0, 1.0e-12));
}

TEST_CASE("IFP ancestry beta_eff accumulator scores delayed ancestor group")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, true);

  const int old_ifp_n_generation = settings::ifp_n_generation;
  const int old_n_inactive = settings::n_inactive;
  settings::ifp_n_generation = 2;
  settings::n_inactive = 10;

  BetaEffectiveAccumulator acc(grid, 2);
  std::unordered_map<int64_t, double> spatial {{0, 1.0}};
  acc.set_adjoint_source_spatial(spatial);

  simulation::clutch_ifp_source_state_bank.clear();
  simulation::clutch_ifp_source_delayed_group_bank.clear();
  simulation::clutch_ifp_source_state_bank.resize(3);
  simulation::clutch_ifp_source_delayed_group_bank.resize(3);
  simulation::clutch_ifp_source_state_bank[0] = {0, 1};
  simulation::clutch_ifp_source_delayed_group_bank[0] = {0, 1};
  simulation::clutch_ifp_source_state_bank[1] = {2, 3};
  simulation::clutch_ifp_source_delayed_group_bank[1] = {2, 0};
  simulation::clutch_ifp_source_state_bank[2] = {4};
  simulation::clutch_ifp_source_delayed_group_bank[2] = {1};

  acc.begin_batch(1);
  acc.score_ifp_ancestry_event(10.0, 0);
  acc.score_ifp_ancestry_event(20.0, 1);
  acc.score_ifp_ancestry_event(30.0, 2);
  acc.end_batch(1);

  auto result = acc.compute_ifp_ancestry_result();
  REQUIRE(result.available);
  REQUIRE_THAT(result.denominator,
    Catch::Matchers::WithinRel(30.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[1],
    Catch::Matchers::WithinRel(20.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[1],
    Catch::Matchers::WithinRel(20.0 / 30.0, 1.0e-12));
  REQUIRE_THAT(result.beta_total,
    Catch::Matchers::WithinRel(20.0 / 30.0, 1.0e-12));

  settings::ifp_n_generation = old_ifp_n_generation;
  settings::n_inactive = old_n_inactive;
}

TEST_CASE("CLUTCH-IFP uses IFP-derived source-state importance")
{
  auto grid = SharedMeshGrid::create(
    1.0, false, {0.0, 0.0, 0.0}, {2.0, 1.0, 1.0}, true);

  const int old_ifp_n_generation = settings::ifp_n_generation;
  const int old_n_inactive = settings::n_inactive;
  settings::ifp_n_generation = 2;
  settings::n_inactive = 10;

  BetaEffectiveAccumulator acc(grid, 2);
  std::unordered_map<int64_t, double> spatial {{0, 1.0}, {9, 1.0}};
  acc.set_adjoint_source_spatial(spatial);

  simulation::clutch_ifp_source_state_bank.clear();
  simulation::clutch_ifp_source_delayed_group_bank.clear();
  simulation::clutch_ifp_source_state_bank.resize(2);
  simulation::clutch_ifp_source_delayed_group_bank.resize(2);
  simulation::clutch_ifp_source_state_bank[0] = {0, 0};
  simulation::clutch_ifp_source_delayed_group_bank[0] = {0, 0};
  simulation::clutch_ifp_source_state_bank[1] = {9, 9};
  simulation::clutch_ifp_source_delayed_group_bank[1] = {0, 0};

  acc.begin_batch(1);
  acc.record_source_birth({0.25, 0.25, 0.25}, 101, 0, 0);
  acc.record_source_birth({1.25, 0.25, 0.25}, 102, 1, 1);
  acc.score_fission_site({0.25, 0.25, 0.25}, 10.0, 0, 1, 1.0, 5.0);
  acc.score_fission_site({1.25, 0.25, 0.25}, 20.0, 1, 2, 1.0, 7.0);
  acc.score_ifp_ancestry_event(4.0, 0);
  acc.score_ifp_ancestry_event(10.0, 1);
  acc.end_batch(1);

  auto result = acc.compute_clutch_ifp_result();
  REQUIRE(result.available);
  REQUIRE_THAT(result.denominator,
    Catch::Matchers::WithinRel(240.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[0],
    Catch::Matchers::WithinRel(40.0, 1.0e-12));
  REQUIRE_THAT(result.numerators[1],
    Catch::Matchers::WithinRel(200.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[0],
    Catch::Matchers::WithinRel(40.0 / 240.0, 1.0e-12));
  REQUIRE_THAT(result.beta_i[1],
    Catch::Matchers::WithinRel(200.0 / 240.0, 1.0e-12));

  settings::ifp_n_generation = old_ifp_n_generation;
  settings::n_inactive = old_n_inactive;
}
