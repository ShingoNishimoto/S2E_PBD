#pragma once

#include "SimulationCase.h"
#include "MCSimExecutor.h"
#include "../Spacecraft/PBD_Sat.h"
#include "RelativeInformation.h"
#include "../InterSatComm/PBD_InterSatComm.h"
#include "../PBD/PBD_pod.h"
#include "../GroundStation/PBD_GroundStation.h"

class PBD_Case : public SimulationCase
{
public:
  PBD_Case(string ini_fname);// , MCSimExecutor& mc_sim, const string log_path);
  virtual ~PBD_Case();
  void Initialize();
  void Main();

  virtual string GetLogHeader() const;
  virtual string GetLogValue() const;

private:
  void InitializeSpacecrafts();
  std::vector<PBD_Sat*> spacecrafts_;
  //MCSimExecutor& mc_sim_;
  RelativeInformation* rel_info_;
  PBD_InterSatComm* pbd_inter_sat_comm_;
  PBD_pod* pbd_;
};
