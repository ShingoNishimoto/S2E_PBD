#pragma once
#include "SimulationConfig.h"

#include "../../Simulation/InterSatComm/RFSystem/RFSystemBeam.h"

class PBD_InterSatComm
{
public:
  PBD_InterSatComm(const SimulationConfig* sim_config);
  ~PBD_InterSatComm();

  // Getter
  RFSystemBeam* GetRFSystemBeam(void) const { return rf_system_beam_; }

private:
  RFSystemBeam* rf_system_beam_;
};

