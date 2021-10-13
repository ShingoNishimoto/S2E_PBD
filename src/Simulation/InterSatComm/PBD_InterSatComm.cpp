#include "../../Interface/InitInput/Initialize_PBD.h"

PBD_InterSatComm::PBD_InterSatComm(const SimulationConfig* sim_config)
{
  rf_system_beam_ = new RFSystemBeam(InitRFSystemBeam(sim_config->inter_sat_comm_file_));
}

PBD_InterSatComm::~PBD_InterSatComm()
{
  delete rf_system_beam_;
}

