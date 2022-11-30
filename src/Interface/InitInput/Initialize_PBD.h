#pragma once
#include <Interface/InitInput/IniAccess.h>
#include <Component/CommGS/InitAntenna.hpp>
#include "../../Component/RFSystem/RFSystemTransmitter.h"
#include "../../Component/RFSystem/RFSystemReceiver.h"
#include "../../Component/AOCS/InitGNSSReceiver.hpp"
#include "../../Component/AOCS/InitializeRelativePositionSensor.hpp"

// Component
class RFSystemTransmitter;
class RFSystemReceiver;
RFSystemTransmitter InitRFSystemTransmitter(ClockGenerator* clock_gen, const std::string ini_path, PBD_InterSatComm* pbd_inter_sat_comm, const Dynamics* dynamics);
RFSystemReceiver InitRFSystemReceiver(ClockGenerator* clock_gen, const std::string ini_path, PBD_InterSatComm* pbd_inter_sat_comm, const Dynamics* dynamics, const SimTime* sim_time);

//InterSatComm
class RFSystemBeam;
RFSystemBeam InitRFSystemBeam(const std::string ini_path);
