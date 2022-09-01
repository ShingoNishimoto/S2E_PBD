#include "PBD_Components.h"
#include "../../Interface/InitInput/Initialize_PBD.h"

PBD_Components::PBD_Components(const Dynamics* dynamics, const Structure* structure, const LocalEnvironment* local_env, const GlobalEnvironment* glo_env, const RelativeInformation* rel_info, PBD_InterSatComm* pbd_inter_sat_comm, const SimulationConfig* config, ClockGenerator* clock_gen, const int sat_id):dynamics_(dynamics), structure_(structure), local_env_(local_env), glo_env_(glo_env), rel_info_(rel_info), pbd_inter_sat_comm_(pbd_inter_sat_comm), config_(config)
{
  IniAccess iniAccess = IniAccess(config->sat_file_[sat_id]);
  // OBC
  obc_ = new OBC(clock_gen);
  // RF System
  std::string ini_path = iniAccess.ReadString("COMPONENTS_FILE", "rf_system_transmitter_file");
  rf_sys_transmitter_ = new RFSystemTransmitter(InitRFSystemTransmitter(clock_gen, ini_path, pbd_inter_sat_comm, dynamics));
  // RF System
  ini_path = iniAccess.ReadString("COMPONENTS_FILE", "rf_system_receiver_file");
  rf_sys_receiver_ = new RFSystemReceiver(InitRFSystemReceiver(clock_gen, ini_path, pbd_inter_sat_comm, dynamics, &(glo_env->GetSimTime())));
  // GNSS receiver
  ini_path = iniAccess.ReadString("COMPONENTS_FILE", "gnss_receiver");
  gnss_receiver_ = new PBD_GNSSReceiver(InitGNSSReceiver(clock_gen, sat_id, ini_path, dynamics, &(glo_env->GetGnssSatellites()), &(glo_env->GetSimTime())));
}

PBD_Components::~PBD_Components()
{
  delete obc_;
  delete rf_sys_transmitter_;
  delete rf_sys_receiver_;
  delete gnss_receiver_;
}

Vector<3> PBD_Components::GenerateForce_N_b()
{
  //There is no orbit control component, so it remains 0
  Vector<3> force_N_b_(0.0);
  return force_N_b_;
};

Vector<3> PBD_Components::GenerateTorque_Nm_b()
{
  //No attitude control component
  Vector<3> torque_Nm_b_(0.0);
  return torque_Nm_b_;
};

void PBD_Components::LogSetUp(Logger & logger)
{
  UNUSED(logger);
}
