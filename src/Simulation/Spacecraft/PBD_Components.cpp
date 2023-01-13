#include "PBD_Components.h"
#include "../../Interface/InitInput/Initialize_PBD.h"

PBD_Components::PBD_Components(const Dynamics* dynamics, const Structure* structure, const LocalEnvironment* local_env, const GlobalEnvironment* glo_env, const RelativeInformation* rel_info, PBD_InterSatComm* pbd_inter_sat_comm, const SimulationConfig* config, ClockGenerator* clock_gen, const int sat_id):dynamics_(dynamics), structure_(structure), local_env_(local_env), glo_env_(glo_env), rel_info_(rel_info), pbd_inter_sat_comm_(pbd_inter_sat_comm), config_(config), sat_id_(sat_id)
{
  IniAccess iniAccess = IniAccess(config->sat_file_[sat_id]);
  double compo_step_sec = glo_env_->GetSimTime().GetCompoStepSec();

  // OBC
  obc_ = new OBC(clock_gen);

  // GNSS receiver
  std::string ini_path = iniAccess.ReadString("COMPONENTS_FILE", "gnss_receiver");
  config_->main_logger_->CopyFileToLogDir(ini_path);
  gnss_receiver_ = new PBD_GNSSReceiver(InitGNSSReceiver(clock_gen, sat_id, ini_path, dynamics, &(glo_env->GetGnssSatellites()), &(glo_env->GetSimTime())));
  PhaseCenterCorrection* pcc_ptr = gnss_receiver_->GetPCCPtr();
  pcc_ptr->LogSetup(*(config_->main_logger_));
  pcc_ptr->PcvLogOutput("pcv_" + std::to_string(sat_id) + ".csv");
  pcc_ptr->PccLogOutput("pcc_" + std::to_string(sat_id) + ".csv");

  if (sat_id == 0)
  {
    // RF System
    ini_path = iniAccess.ReadString("COMPONENTS_FILE", "rf_system_transmitter_file");
    rf_sys_transmitter_ = new RFSystemTransmitter(InitRFSystemTransmitter(clock_gen, ini_path, pbd_inter_sat_comm, dynamics));

    // RF System
    ini_path = iniAccess.ReadString("COMPONENTS_FILE", "rf_system_receiver_file");
    rf_sys_receiver_ = new RFSystemReceiver(InitRFSystemReceiver(clock_gen, ini_path, pbd_inter_sat_comm, dynamics, &(glo_env->GetSimTime())));

    // Relative Position Sensor
    const std::string rel_pos_file = iniAccess.ReadString("COMPONENTS_FILE", "relative_position_sensor_file");
    config_->main_logger_->CopyFileToLogDir(rel_pos_file);
    relative_position_sensor_ = new RelativePositionSensor(InitializeRelativePositionSensor(clock_gen, rel_pos_file, compo_step_sec, *rel_info_, *dynamics_));
  }
}

PBD_Components::~PBD_Components()
{
  if (sat_id_ == 0)
  {
    delete rf_sys_transmitter_;
    delete rf_sys_receiver_;
    delete relative_position_sensor_;
  }

  delete gnss_receiver_;
  // OBC must be deleted the last since it has com ports
  delete obc_;
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

void PBD_Components::LogSetup(Logger & logger)
{
  logger.AddLoggable(gnss_receiver_);
  if (sat_id_ == 0)
  {
    logger.AddLoggable(relative_position_sensor_);
  }
}
