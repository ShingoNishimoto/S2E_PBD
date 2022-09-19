#pragma once

#include "PBD_GNSSReceiver.hpp"

PBD_GNSSReceiver InitGNSSReceiver(ClockGenerator* clock_gen, int id, const std::string fname,
  const Dynamics* dynamics, const GnssSatellites* gnss_satellites, const SimTime* simtime);
