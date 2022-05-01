#pragma once

#include "Spacecraft.h"
#include "PBD_Components.h"
#include "../PBD/PBD_GnssObservation.h"

class PBD_Components;

class PBD_Sat : public Spacecraft
{
public:
  PBD_Sat(SimulationConfig* sim_config, const GlobalEnvironment* glo_env, RelativeInformation* rel_info, PBD_InterSatComm* pbd_inter_sat_comm, const int sat_id);
  ~PBD_Sat();

  // 初期化
  virtual void Initialize(SimulationConfig* sim_config, const GlobalEnvironment* glo_env, PBD_InterSatComm* pbd_inter_sat_comm, const int sat_id);
  // ログ保存機能
  virtual void LogSetup(Logger& logger);
  // 状態量の更新
  virtual void Update(const SimTime* sim_time);

  //ダイナミクスへの力・トルク出力
  void GenerateTorque_b();
  void GenerateForce_b();

  PBD_GnssObservation& gnss_observation_;

private:
  PBD_Components* components_;
};