#pragma once

#include <SimulationCase.h>
#include <MCSimExecutor.h>
#include <RelativeInformation.h>
#include "../Spacecraft/PBD_Sat.h"
#include "../InterSatComm/PBD_InterSatComm.h"
#include "../PBD/PBD_dgps.h"
#include "../GroundStation/PBD_GroundStation.h"

class PBD_Case : public SimulationCase
{
public:
  PBD_Case(std::string ini_base);// , MCSimExecutor& mc_sim, const string log_path);
  virtual ~PBD_Case();

  void Initialize();
  void Main();

  virtual std::string GetLogHeader() const;
  virtual std::string GetLogValue() const;

private:
  void InitializeSpacecrafts();
  std::vector<PBD_Sat*> spacecrafts_;
  //MCSimExecutor& mc_sim_;
  RelativeInformation rel_info_;
  PBD_InterSatComm* pbd_inter_sat_comm_;
  PBD_dgps* pbd_;
};
