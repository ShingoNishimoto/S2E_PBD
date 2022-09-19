#include "PBD_Sat.h"

PBD_Sat::PBD_Sat(SimulationConfig* sim_config, const GlobalEnvironment* glo_env, RelativeInformation* rel_info, PBD_InterSatComm* pbd_inter_sat_comm, const int sat_id)
:Spacecraft(sim_config, glo_env, rel_info, sat_id)//, gnss_observation_(PBD_GnssObservation(this->GetDynamics().GetOrbit(), glo_env->GetGnssSatellites()))
{
  // Initialize(sim_config, glo_env, pbd_inter_sat_comm, sat_id);
  pbd_components_ = new PBD_Components(dynamics_, structure_, local_env_, glo_env, rel_info_, pbd_inter_sat_comm, sim_config, &clock_gen_, sat_id);
  components_ = pbd_components_; // これでいいのか？
  gnss_observation_ = new PBD_GnssObservation(pbd_components_->GetGNSSReceiver(), glo_env->GetGnssSatellites());
}

PBD_Sat::~PBD_Sat()
{
  delete gnss_observation_;
}

/*
void PBD_Sat::Initialize(SimulationConfig* sim_config, const GlobalEnvironment* glo_env, PBD_InterSatComm* pbd_inter_sat_comm, const int sat_id)
{
  components_ = new PBD_Components(dynamics_, structure_, local_env_, glo_env, rel_info_, pbd_inter_sat_comm, sim_config, &clock_gen_, sat_id);
}
*/

void PBD_Sat::LogSetup(Logger & logger)
{
  Spacecraft::LogSetup(logger);
  // pbd_components_->PBD_Components::LogSetup(logger);
}

void PBD_Sat::Update(const SimTime* sim_time)
{
  Spacecraft::Update(sim_time);
  // この更新を毎回すると計算は重そうだが，しょうがないか．<- 実際の観測頻度に合わせてもいいかも
  gnss_observation_->Update();
}
