#pragma once

#include "GroundStation.h"
#include "PBD_GSComponents.h"
#include "../Spacecraft/PBD_Sat.h"


// class PBD_GSComponents;

class PBD_GroundStation : public GroundStation
{
public:
  PBD_GroundStation(SimulationConfig* config, int gs_id);
  ~PBD_GroundStation();

  // 初期化
  virtual void Initialize(SimulationConfig* config);
  // ログ保存機能
  virtual void LogSetup(Logger& logger);
  // 状態量の更新
  virtual void Update(const Dynamics& dynamics, const ANT& sc_ant, const PBD_GroundStation& PBDgs, const double& current_jd);

private:
  PBD_GSComponents* components_;
};
