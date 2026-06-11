#ifndef OPENMC_CLUTCH_IFP_H
#define OPENMC_CLUTCH_IFP_H

#include "openmc/ifp.h"
#include "openmc/message_passing.h"
#include "openmc/particle.h"
#include "openmc/vector.h"

#include <cstdint>

namespace openmc {

bool clutch_ifp_on();

int clutch_ifp_n_generation();

void resize_simulation_clutch_ifp_banks();

void record_clutch_ifp_fission_site(
  const Particle& p, int64_t child_state, int delayed_group, int64_t idx);

void copy_clutch_ifp_data_from_fission_banks(int i_bank,
  vector<int64_t>& states, vector<int>& delayed_groups);

void allocate_temporary_vector_clutch_ifp(
  vector<vector<int64_t>>& states, vector<vector<int>>& delayed_groups);

void copy_clutch_ifp_data_to_fission_banks(
  const vector<int64_t>* states_ptr, const vector<int>* delayed_groups_ptr);

#ifdef OPENMC_MPI

void broadcast_clutch_ifp_n_generation(int& n_generation,
  const vector<vector<int64_t>>& states,
  const vector<vector<int>>& delayed_groups);

void send_clutch_ifp_info(int64_t idx, int64_t n, int n_generation,
  int neighbor, vector<MPI_Request>& requests,
  const vector<vector<int64_t>>& states, vector<int64_t>& send_states,
  const vector<vector<int>>& delayed_groups, vector<int>& send_delayed_groups);

void receive_clutch_ifp_data(int64_t idx, int64_t n, int n_generation,
  int neighbor, vector<MPI_Request>& requests, vector<int64_t>& states,
  vector<int>& delayed_groups, vector<DeserializationInfo>& deserialization);

void copy_partial_clutch_ifp_data_to_source_banks(int64_t idx, int n,
  int64_t i_bank, const vector<vector<int64_t>>& states,
  const vector<vector<int>>& delayed_groups);

void deserialize_clutch_ifp_info(int n_generation,
  const vector<DeserializationInfo>& deserialization,
  const vector<int64_t>& states, const vector<int>& delayed_groups);

#endif

void copy_complete_clutch_ifp_data_to_source_banks(
  const vector<vector<int64_t>>& states,
  const vector<vector<int>>& delayed_groups);

} // namespace openmc

#endif // OPENMC_CLUTCH_IFP_H
