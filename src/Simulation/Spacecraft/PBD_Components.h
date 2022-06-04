#pragma once
#include "Vector.hpp"
#include "Dynamics.h"
#include "RelativeInformation.h"
#include "Structure.h"
#include "../../Simulation/InterSatComm/PBD_InterSatComm.h"
#include <GlobalEnvironment.h>
#include <LocalEnvironment.h>
#include <Simulation/Spacecraft/InstalledComponents.hpp>

#include "../../Component/CDH/OBC_Sat0.h"
#include "../../Component/CDH/OBC_Sat1.h"
#include "../../Component/RFSystem/RFSystemTransmitter.h"
#include "../../Component/RFSystem/RFSystemReceiver.h"

//class OBC_Sat0;

using libra::Vector;

class PBD_Components : public InstalledComponents
{
public:
  PBD_Components(
    const Dynamics* dynamics,
    const Structure* structure,
    const LocalEnvironment* local_env,
    const GlobalEnvironment* glo_env,
    const RelativeInformation* rel_info,
    PBD_InterSatComm* pbd_inter_sat_comm,
    const SimulationConfig* config,
    ClockGenerator* clock_gen,
    const int sat_id);
  ~PBD_Components();
  Vector<3> GenerateForce_N_b();
  Vector<3> GenerateTorque_Nm_b();
  void LogSetUp(Logger& logger);
  //TODO: Do null-check in the getter (Avoid pointing to another satellite component)

private:
  // この構成はなんか変な気がするな．．
  // Sat0
  OBC_Sat0* obc0_;
  RFSystemTransmitter* rf_sys_transmitter_;

  // Sat1
  OBC_Sat1* obc1_;
  RFSystemReceiver* rf_sys_receiver_;

  // General information
  const SimulationConfig* config_;
  const Dynamics* dynamics_;
  const Structure* structure_;
  const LocalEnvironment* local_env_;
  const GlobalEnvironment* glo_env_;
  const RelativeInformation* rel_info_;
  PBD_InterSatComm* pbd_inter_sat_comm_;
  const int sat_id_;
};
