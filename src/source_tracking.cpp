#include "openmc/source_tracking.h"

#include <fstream>
#include <iostream>
#include <vector>

#include "openmc/particle.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

namespace openmc {

//==============================================================================
// Static variables
//==============================================================================

static bool tracking_enabled_ = false;
static std::vector<SourceTrackingInfo> tracking_data_;
static std::ofstream tracking_file_;

//==============================================================================
// Functions
//==============================================================================

void init_source_tracking()
{
  // Enable tracking only for active generations in eigenvalue mode
  if (settings::run_mode == RunMode::EIGENVALUE &&
      simulation::current_batch > settings::n_inactive) {
    tracking_enabled_ = true;

    // Open output file
    std::string filename = "source_tracking_batch_" +
                           std::to_string(simulation::current_batch) + ".txt";
    tracking_file_.open(filename);

    if (tracking_file_.is_open()) {
      // Write header
      tracking_file_ << "# Source Particle Tracking Data\n";
      tracking_file_
        << "# Columns: source_label source_batch source_x source_y source_z "
        << "current_x current_y current_z energy weight collisions is_source\n";
    }

    tracking_data_.clear();
  } else if (settings::run_mode == RunMode::FIXED_SOURCE) {
    tracking_enabled_ = true;

    std::string filename =
      "source_tracking_gen_" + std::to_string(simulation::total_gen) + ".txt";
    tracking_file_.open(filename);

    if (tracking_file_.is_open()) {
      tracking_file_ << "# Source Particle Tracking Data (Fixed Source)\n";
      tracking_file_
        << "# Columns: source_label generation source_x source_y source_z "
        << "current_x current_y current_z energy weight collisions is_source\n";
    }

    tracking_data_.clear();
  }
}

void record_source_tracking(const Particle& p)
{
  if (!tracking_enabled_ || !tracking_file_.is_open())
    return;

  // Only record particles with valid source labels
  if (p.source_label() == 0)
    return;

  SourceTrackingInfo info;
  info.source_label = p.source_label();
  info.source_position = p.source_position();
  info.source_batch = p.source_batch();
  info.current_position = p.r();
  info.current_energy = p.E();
  info.weight = p.wgt();
  info.collision_count = p.n_collision();
  info.is_source_particle =
    ((p.r() - p.source_position()).norm() < 1e-6 && p.n_collision() == 0);

  // Write to file immediately
  tracking_file_ << info.source_label << " " << info.source_batch << " "
                 << info.source_position.x << " " << info.source_position.y
                 << " " << info.source_position.z << " "
                 << info.current_position.x << " " << info.current_position.y
                 << " " << info.current_position.z << " " << info.current_energy
                 << " " << info.weight << " " << info.collision_count << " "
                 << (info.is_source_particle ? 1 : 0) << "\n";
}

void finalize_source_tracking()
{
  if (tracking_file_.is_open()) {
    tracking_file_.close();
  }
  tracking_enabled_ = false;
  tracking_data_.clear();
}

bool source_tracking_enabled()
{
  return tracking_enabled_;
}

} // namespace openmc
