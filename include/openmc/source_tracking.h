#ifndef OPENMC_SOURCE_TRACKING_H
#define OPENMC_SOURCE_TRACKING_H

#include "openmc/particle_data.h"
#include "openmc/position.h"

namespace openmc {

// Forward declarations
class Particle;

//! Structure to hold source particle tracking information
struct SourceTrackingInfo {
  int64_t source_label;      //!< Unique source particle label
  Position source_position;  //!< Initial source position
  int source_batch;          //!< Batch when source was created
  Position current_position; //!< Current particle position
  double current_energy;     //!< Current particle energy
  double weight;             //!< Current particle weight
  int collision_count;       //!< Number of collisions
  bool is_source_particle;   //!< True if this is the original source particle
};

//! Initialize source particle tracking system
void init_source_tracking();

//! Record source particle tracking information
void record_source_tracking(const Particle& p);

//! Finalize and write source tracking data
void finalize_source_tracking();

//! Check if source tracking is enabled
bool source_tracking_enabled();

} // namespace openmc

#endif // OPENMC_SOURCE_TRACKING_H
