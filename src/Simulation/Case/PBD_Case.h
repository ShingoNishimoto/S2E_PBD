#pragma once

#include "SimulationCase.h"
#include "MCSimExecutor.h"
#include "../Spacecraft/PBD_Sat.h"
#include "RelativeInformation.h"
#include "../InterSatComm/PBD_InterSatComm.h"

class PBD_Case : public SimulationCase
{
public:
  PBD_Case(string ini_fname);
  virtual ~PBD_Case();
  
  void Initialize();
  void Main();

  virtual string GetLogHeader() const;
  virtual string GetLogValue() const;

private:
  void InitializeSpacecrafts();
  std::vector<PBD_Sat*> spacecrafts_;
  RelativeInformation* rel_info_;
  PBD_InterSatComm* pbd_inter_sat_comm_;
};
