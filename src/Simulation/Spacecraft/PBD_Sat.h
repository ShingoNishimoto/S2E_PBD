#pragma once

#include <Spacecraft.h>
#include "PBD_Components.h"
#include "../PBD/PBD_GnssObservation.h"

class PBD_Sat : public Spacecraft
{
public:
  PBD_Sat(SimulationConfig* sim_config, const GlobalEnvironment* glo_env, RelativeInformation* rel_info, PBD_InterSatComm* pbd_inter_sat_comm, const int sat_id);

  // 状態量の更新
  virtual void Update(const SimTime* sim_time);
  ~PBD_Sat();

  // ログ保存機能
  virtual void LogSetup(Logger& logger);


  // FIXME: この辺は直したい．
  PBD_GnssObservation* gnss_observation_;
  PBD_Components* pbd_components_; // これはprivateにすべきかも．

private:
};