#include "openmc/clutch_ifp.h"

#include "openmc/bank.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

#include <algorithm>

namespace openmc {

namespace {

template<typename T>
vector<T> push_clutch_ifp_value(
  const T& value, const vector<T>& data, int n_generation)
{
  vector<T> updated;
  const size_t source_idx = data.size();

  if (source_idx < static_cast<size_t>(n_generation)) {
    updated.resize(source_idx + 1);
    for (size_t i = 0; i < source_idx; ++i) {
      updated[i] = data[i];
    }
    updated[source_idx] = value;
  } else if (source_idx == static_cast<size_t>(n_generation)) {
    updated.resize(source_idx);
    for (size_t i = 0; i < source_idx - 1; ++i) {
      updated[i] = data[i + 1];
    }
    updated[source_idx - 1] = value;
  }
  return updated;
}

} // namespace

bool clutch_ifp_on()
{
  return settings::beta_effective_on &&
         settings::run_mode == RunMode::EIGENVALUE &&
         settings::solver_type == SolverType::MONTE_CARLO;
}

int clutch_ifp_n_generation()
{
  int n_generation = settings::ifp_n_generation > 0 ?
                       settings::ifp_n_generation :
                       DEFAULT_IFP_N_GENERATION;
  if (settings::n_inactive > 0) {
    n_generation = std::min(n_generation, settings::n_inactive);
  }
  return std::max(1, n_generation);
}

void resize_simulation_clutch_ifp_banks()
{
  if (!clutch_ifp_on()) {
    return;
  }
  simulation::clutch_ifp_source_state_bank.resize(simulation::work_per_rank);
  simulation::clutch_ifp_source_delayed_group_bank.resize(
    simulation::work_per_rank);
  simulation::clutch_ifp_fission_state_bank.resize(
    3 * simulation::work_per_rank);
  simulation::clutch_ifp_fission_delayed_group_bank.resize(
    3 * simulation::work_per_rank);
}

void record_clutch_ifp_fission_site(
  const Particle& p, int64_t child_state, int delayed_group, int64_t idx)
{
  if (!clutch_ifp_on() || idx < 0) {
    return;
  }

  const int64_t source_index = p.current_work() - 1;
  const int n_generation = clutch_ifp_n_generation();

  vector<int64_t> parent_states;
  vector<int> parent_delayed_groups;
  if (source_index >= 0 &&
      source_index <
        static_cast<int64_t>(simulation::clutch_ifp_source_state_bank.size())) {
    parent_states = simulation::clutch_ifp_source_state_bank[source_index];
    parent_delayed_groups =
      simulation::clutch_ifp_source_delayed_group_bank[source_index];
  }

  simulation::clutch_ifp_fission_state_bank[idx] =
    push_clutch_ifp_value(child_state, parent_states, n_generation);
  simulation::clutch_ifp_fission_delayed_group_bank[idx] =
    push_clutch_ifp_value(delayed_group, parent_delayed_groups, n_generation);
}

void copy_clutch_ifp_data_from_fission_banks(int i_bank,
  vector<int64_t>& states, vector<int>& delayed_groups)
{
  if (!clutch_ifp_on()) {
    return;
  }
  states = simulation::clutch_ifp_fission_state_bank[i_bank];
  delayed_groups = simulation::clutch_ifp_fission_delayed_group_bank[i_bank];
}

void allocate_temporary_vector_clutch_ifp(
  vector<vector<int64_t>>& states, vector<vector<int>>& delayed_groups)
{
  if (!clutch_ifp_on()) {
    return;
  }
  states.resize(simulation::fission_bank.size());
  delayed_groups.resize(simulation::fission_bank.size());
}

void copy_clutch_ifp_data_to_fission_banks(
  const vector<int64_t>* states_ptr, const vector<int>* delayed_groups_ptr)
{
  if (!clutch_ifp_on()) {
    return;
  }
  std::copy(states_ptr, states_ptr + simulation::fission_bank.size(),
    simulation::clutch_ifp_fission_state_bank.data());
  std::copy(delayed_groups_ptr,
    delayed_groups_ptr + simulation::fission_bank.size(),
    simulation::clutch_ifp_fission_delayed_group_bank.data());
}

#ifdef OPENMC_MPI

void broadcast_clutch_ifp_n_generation(int& n_generation,
  const vector<vector<int64_t>>& states,
  const vector<vector<int>>& delayed_groups)
{
  if (mpi::rank == 0) {
    if (!states.empty()) {
      n_generation = static_cast<int>(states[0].size());
    } else if (!delayed_groups.empty()) {
      n_generation = static_cast<int>(delayed_groups[0].size());
    } else {
      n_generation = clutch_ifp_n_generation();
    }
  }
  MPI_Bcast(&n_generation, 1, MPI_INT, 0, mpi::intracomm);
}

void send_clutch_ifp_info(int64_t idx, int64_t n, int n_generation,
  int neighbor, vector<MPI_Request>& requests,
  const vector<vector<int64_t>>& states, vector<int64_t>& send_states,
  const vector<vector<int>>& delayed_groups, vector<int>& send_delayed_groups)
{
  for (int64_t i = idx; i < idx + n; ++i) {
    std::copy(states[i].begin(), states[i].end(),
      send_states.begin() + i * n_generation);
    std::copy(delayed_groups[i].begin(), delayed_groups[i].end(),
      send_delayed_groups.begin() + i * n_generation);
  }

  requests.emplace_back();
  MPI_Isend(&send_states[n_generation * idx],
    n_generation * static_cast<int>(n), MPI_INT64_T, neighbor, mpi::rank,
    mpi::intracomm, &requests.back());

  requests.emplace_back();
  MPI_Isend(&send_delayed_groups[n_generation * idx],
    n_generation * static_cast<int>(n), MPI_INT, neighbor, mpi::rank,
    mpi::intracomm, &requests.back());
}

void receive_clutch_ifp_data(int64_t idx, int64_t n, int n_generation,
  int neighbor, vector<MPI_Request>& requests, vector<int64_t>& states,
  vector<int>& delayed_groups, vector<DeserializationInfo>& deserialization)
{
  requests.emplace_back();
  MPI_Irecv(&states[n_generation * idx], n_generation * static_cast<int>(n),
    MPI_INT64_T, neighbor, neighbor, mpi::intracomm, &requests.back());

  requests.emplace_back();
  MPI_Irecv(&delayed_groups[n_generation * idx],
    n_generation * static_cast<int>(n), MPI_INT, neighbor, neighbor,
    mpi::intracomm, &requests.back());

  DeserializationInfo info = {idx, n};
  deserialization.push_back(info);
}

void copy_partial_clutch_ifp_data_to_source_banks(int64_t idx, int n,
  int64_t i_bank, const vector<vector<int64_t>>& states,
  const vector<vector<int>>& delayed_groups)
{
  std::copy(&states[idx], &states[idx + n],
    &simulation::clutch_ifp_source_state_bank[i_bank]);
  std::copy(&delayed_groups[idx], &delayed_groups[idx + n],
    &simulation::clutch_ifp_source_delayed_group_bank[i_bank]);
}

void deserialize_clutch_ifp_info(int n_generation,
  const vector<DeserializationInfo>& deserialization,
  const vector<int64_t>& states, const vector<int>& delayed_groups)
{
  for (auto info : deserialization) {
    int64_t index_local = info.index_local;
    int64_t n = info.n;

    for (int64_t i = index_local; i < index_local + n; ++i) {
      vector<int64_t> states_received(
        states.begin() + n_generation * i,
        states.begin() + n_generation * (i + 1));
      vector<int> delayed_groups_received(
        delayed_groups.begin() + n_generation * i,
        delayed_groups.begin() + n_generation * (i + 1));
      simulation::clutch_ifp_source_state_bank[i] = states_received;
      simulation::clutch_ifp_source_delayed_group_bank[i] =
        delayed_groups_received;
    }
  }
}

#endif

void copy_complete_clutch_ifp_data_to_source_banks(
  const vector<vector<int64_t>>& states,
  const vector<vector<int>>& delayed_groups)
{
  if (!clutch_ifp_on()) {
    return;
  }
  std::copy(states.data(), states.data() + settings::n_particles,
    simulation::clutch_ifp_source_state_bank.begin());
  std::copy(delayed_groups.data(), delayed_groups.data() + settings::n_particles,
    simulation::clutch_ifp_source_delayed_group_bank.begin());
}

} // namespace openmc
