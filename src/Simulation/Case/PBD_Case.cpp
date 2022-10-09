#include "PBD_Case.h"
#include <Interface/InitInput/IniAccess.h>
#include "../PBD/PBD_GeoPotential.hpp"

PBD_Case::PBD_Case(std::string ini_base) :SimulationCase(ini_base)//, MCSimExecutor& mc_sim, const string log_path):SimulationCase(ini_fname, mc_sim, log_path),mc_sim_(mc_sim)
{
  rel_info_ = RelativeInformation();
  pbd_inter_sat_comm_ = new PBD_InterSatComm(&sim_config_);
}

PBD_Case::~PBD_Case()
{
  for (auto& spacefraft : spacecrafts_)
  {
    delete spacefraft;
  }
  delete pbd_inter_sat_comm_;
}

// Spacecraft initializer that manage the order in which each satellite is initialized
// NOTE: Satellites using RelativeOrbit need to be initialized after the reference satellite initialization.
void PBD_Case::InitializeSpacecrafts()
{
  std::vector<int> relative_orbit_sat_id; // List of satellite IDs that use RelativeOrbit

  for (int sat_id = 0; sat_id < sim_config_.num_of_simulated_spacecraft_; sat_id++)
  {
    PBD_Sat* spacecraft = new PBD_Sat(&sim_config_, glo_env_, &rel_info_, pbd_inter_sat_comm_, sat_id);
    spacecrafts_.push_back(spacecraft);
  }
}

void PBD_Case::Initialize()
{
  InitializeSpacecrafts();

  //Register the log output
  glo_env_->LogSetup(*(sim_config_.main_logger_));
  for (auto& spacecraft : spacecrafts_)
  {
    spacecraft->LogSetup(*(sim_config_.main_logger_));
  }

  // この時点で受け取っても中身がないのかもしれない．なぜか中身がなくて怒られてしまう．
  const LocalEnvironment& local_env_ = spacecrafts_.at(0)->GetLocalEnv();
  PBD_GeoPotential geop(2, local_env_, "../../../ExtLibraries/GeoPotential/egm96_to360.ascii");

  pbd_ = new PBD_dgps(glo_env_->GetSimTime(), glo_env_->GetGnssSatellites(), spacecrafts_.at(0)->GetDynamics(), spacecrafts_.at(1)->GetDynamics(), *(spacecrafts_.at(0)->gnss_observation_), *(spacecrafts_.at(1)->gnss_observation_), geop); // ここはGetterとか使った方がいい．

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
      // コンポーネントの更新もここでやるのか？MainRoutineがどこで呼ばれているのかが不明．．．

    }
    // ここで毎回local Environmentを受け取ったほうがいいのか？
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

  return str_tmp;
}

std::string PBD_Case::GetLogValue() const
{
  std::string str_tmp = "";

  return str_tmp;
}
