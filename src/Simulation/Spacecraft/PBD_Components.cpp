#include "PBD_Components.h"
#include "../../Interface/InitInput/Initialize_PBD.h"

PBD_Components::PBD_Components(const Dynamics* dynamics, const Structure* structure, const LocalEnvironment* local_env, const GlobalEnvironment* glo_env, const RelativeInformation* rel_info, PBD_InterSatComm* pbd_inter_sat_comm, const SimulationConfig* config, ClockGenerator* clock_gen, const int sat_id):dynamics_(dynamics), structure_(structure), local_env_(local_env), glo_env_(glo_env), rel_info_(rel_info), pbd_inter_sat_comm_(pbd_inter_sat_comm), config_(config), sat_id_(sat_id)
{
  IniAccess iniAccess = IniAccess(config->sat_file_[sat_id]);
  if(sat_id_ == 0)
  {
    // OBC
    obc0_ = new OBC_Sat0(clock_gen, *this);
    // RF System
    string ini_path = iniAccess.ReadString("COMPONENTS_FILE", "rf_system_transmitter_file");
    rf_sys_transmitter_ = new RFSystemTransmitter(InitRFSystemTransmitter(clock_gen, ini_path, pbd_inter_sat_comm, dynamics));
  }
  else if (sat_id_ == 1)
  {
    // OBC
    obc1_ = new OBC_Sat1(clock_gen, *this);
    // RF System
    string ini_path = iniAccess.ReadString("COMPONENTS_FILE", "rf_system_receiver_file");
    rf_sys_receiver_ = new RFSystemReceiver(InitRFSystemReceiver(clock_gen, ini_path, pbd_inter_sat_comm, dynamics, &(glo_env->GetSimTime())));
  }
  else
  {
    // not reached
  }
}

PBD_Components::~PBD_Components()
{
  switch (sat_id_)
  {
  case 0:
    delete obc0_;
    delete rf_sys_transmitter_;
    break;
  case 1:
    delete obc1_;
    delete rf_sys_receiver_;
    break;
  }
}

Vector<3> PBD_Components::GenerateForce_b()
{
  //There is no orbit control component, so it remains 0
  Vector<3> force_b_(0.0);
  return force_b_;
};

Vector<3> PBD_Components::GenerateTorque_b()
{
  //No attitude control component
  Vector<3> torque_b_(0.0);
  return torque_b_;
};

void PBD_Components::CompoLogSetUp(Logger & logger)
{
}