#pragma once

#include <Simulation/Spacecraft/InstalledComponents.hpp>

#include "Vector.hpp"
#include "Dynamics.h"
#include <GlobalEnvironment.h>
#include <LocalEnvironment.h>

// include for components
#include <Component/CDH/OBC.h>
#include "../../Component/RFSystem/RFSystemTransmitter.h"
#include "../../Component/RFSystem/RFSystemReceiver.h"
#include "../../Component/AOCS/PBD_GNSSReceiver.hpp"
#include "../../Component/AOCS/RelativePositionSensor.hpp"

class PBD_Components : public InstalledComponents
{
public:
  PBD_Components(const Dynamics* dynamics, const Structure* structure, const LocalEnvironment* local_env,
  const GlobalEnvironment* glo_env, const RelativeInformation* rel_info, PBD_InterSatComm* pbd_inter_sat_comm,
  const SimulationConfig* config, ClockGenerator* clock_gen, const int sat_id);
  ~PBD_Components();
  libra::Vector<3> GenerateForce_N_b();
  libra::Vector<3> GenerateTorque_Nm_b();
  void LogSetup(Logger& logger) override;

  // Getter
  //TODO: Do null-check in the getter (Avoid pointing to another satellite component)
  PBD_GNSSReceiver* GetGNSSReceiver(void) const { return gnss_receiver_; }
  RelativePositionSensor* GetRelativePositionSensor(void) const { return relative_position_sensor_; }

private:
  const int sat_id_;

  OBC* obc_;
  RFSystemTransmitter* rf_sys_transmitter_;
  RFSystemReceiver* rf_sys_receiver_;
  PBD_GNSSReceiver* gnss_receiver_;
  RelativePositionSensor* relative_position_sensor_;

  // General information
  const SimulationConfig* config_;
  const Dynamics* dynamics_;
  const Structure* structure_;
  const LocalEnvironment* local_env_;
  const GlobalEnvironment* glo_env_;
  const RelativeInformation* rel_info_;
  PBD_InterSatComm* pbd_inter_sat_comm_;
};
