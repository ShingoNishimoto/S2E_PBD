#pragma once
#include "Initialize.h"
#include "../../Simulation/InterSatComm/PBD_InterSatComm.h"

// Component
class RFSystemTransmitter;
class RFSystemReceiver;
RFSystemTransmitter InitRFSystemTransmitter(ClockGenerator* clock_gen, const string ini_path, PBD_InterSatComm* pbd_inter_sat_comm, const Dynamics* dynamics);
RFSystemReceiver InitRFSystemReceiver(ClockGenerator* clock_gen, const string ini_path, PBD_InterSatComm* pbd_inter_sat_comm, const Dynamics* dynamics, const SimTime* sim_time);

//InterSatComm
class RFSystemBeam;
RFSystemBeam InitRFSystemBeam(const string ini_path);
