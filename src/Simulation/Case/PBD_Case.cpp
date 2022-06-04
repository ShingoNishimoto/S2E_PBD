#include "PBD_Case.h"
#include <Interface/InitInput/IniAccess.h>

PBD_Case::PBD_Case(std::string ini_fname) :SimulationCase(ini_fname)//, MCSimExecutor& mc_sim, const string log_path):SimulationCase(ini_fname, mc_sim, log_path),mc_sim_(mc_sim)
{
  rel_info_ = new RelativeInformation();
  pbd_inter_sat_comm_ = new PBD_InterSatComm(&sim_config_);
}

PBD_Case::~PBD_Case()
{
  for (auto& spacefraft : spacecrafts_)
  {
    delete spacefraft;
  }
  delete rel_info_;
  delete pbd_inter_sat_comm_;
}

// Spacecraft initializer that manage the order in which each satellite is initialized
// NOTE: Satellites using RelativeOrbit need to be initialized after the reference satellite initialization.
void PBD_Case::InitializeSpacecrafts()
{
  std::vector<int> relative_orbit_sat_id; // List of satellite IDs that use RelativeOrbit

  for (int sat_id = 0; sat_id < sim_config_.num_of_simulated_spacecraft_; sat_id++)
  {
    auto orbit_conf = IniAccess(sim_config_.sat_file_[sat_id]); // FIXME: ここにあるのは使いにくい．
    char* section = "ORBIT";
    Orbit::PROPAGATE_MODE propagate_mode = Orbit::PROPAGATE_MODE(orbit_conf.ReadInt(section, "propagate_mode"));

    if (propagate_mode == Orbit::PROPAGATE_MODE::RELATIVE_ORBIT)
    {
      // Memorize the IDs of satellites that use RelativeOrbit, but do not instantiate FFSat
      relative_orbit_sat_id.push_back(sat_id);
    }
    else
    {
      // For satellites that do not use RelativeOrbit, instantiate as usual
      PBD_Sat* spacecraft = new PBD_Sat(&sim_config_, glo_env_, rel_info_, pbd_inter_sat_comm_, sat_id);
      spacecrafts_.push_back(spacecraft);
    }
  }

  // Then, instantiate satellites that use RelativeOrbit
  for (auto sat_id : relative_orbit_sat_id)
  {
    PBD_Sat* spacecraft = new PBD_Sat(&sim_config_, glo_env_, rel_info_, pbd_inter_sat_comm_, sat_id);
    spacecrafts_.push_back(spacecraft);
  }
}

void PBD_Case::Initialize()
{
  InitializeSpacecrafts();
  pbd_ = new PBD_dgps(glo_env_->GetSimTime(), glo_env_->GetGnssSatellites(), (spacecrafts_.at(0)->GetDynamics()).GetOrbit(), (spacecrafts_.at(1)->GetDynamics()).GetOrbit(), *(spacecrafts_.at(0)->gnss_observation_), *(spacecrafts_.at(1)->gnss_observation_)); // ここはGetterとか使った方がいい．
  //Register the log output
  glo_env_->LogSetup(*(sim_config_.main_logger_));
  rel_info_->LogSetup(*(sim_config_.main_logger_));
  for (auto& spacecraft : spacecrafts_)
  {
    spacecraft->LogSetup(*(sim_config_.main_logger_));
  }

  //Write headers to the log
  sim_config_.main_logger_->WriteHeaders();

  //Start the simulation
  std::cout << "\nSimulationDateTime \n";
  glo_env_->GetSimTime().PrintStartDateTime();

}

void PBD_Case::Main()
{
  glo_env_->Reset(); //for MonteCarlo Sim
  while (!glo_env_->GetSimTime().GetState().finish)
  {
    //Logging
    if (glo_env_->GetSimTime().GetState().log_output)
    {
      sim_config_.main_logger_->WriteValues();
    }
    // Global Environment Update
    glo_env_->Update();
    // Spacecraft Update
    for (auto& spacecraft : spacecrafts_)
    {
      // 実際はここの中でGNSS観測情報の更新をしたい．
      spacecraft->Update(&(glo_env_->GetSimTime()));
      spacecraft->Clear(); //Zero clear force and torque for dynamics
    }
    // 軌道のupdateができてるかどうか確認した方がいい
    pbd_->Update(glo_env_->GetSimTime(), glo_env_->GetGnssSatellites(), *(spacecrafts_.at(0)->gnss_observation_), *(spacecrafts_.at(1)->gnss_observation_));
    // Debug output
    if (glo_env_->GetSimTime().GetState().disp_output)
    {
      std::cout << "Progresss: " << glo_env_->GetSimTime().GetProgressionRate() << "%\r";
    }
  }
}

std::string PBD_Case::GetLogHeader() const
{
  std::string str_tmp = "";
  str_tmp += WriteScalar("time", "s");
  for (auto& spacecraft : spacecrafts_)
  {
    str_tmp += WriteVector("Sat" + std::to_string(spacecraft->GetSatID()) + "_Omega", "b", "rad/s", 3);
  }
  return str_tmp;
}

std::string PBD_Case::GetLogValue() const
{
  std::string str_tmp = "";
  str_tmp += WriteScalar(glo_env_->GetSimTime().GetElapsedSec());
  for (auto& spacecraft : spacecrafts_)
  {
    str_tmp += WriteVector(spacecraft->GetDynamics().GetAttitude().GetOmega_b(), 3);
  }
  return str_tmp;
}
