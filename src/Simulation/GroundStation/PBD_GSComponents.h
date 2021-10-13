#pragma once
#include "ANT.h"
#include "GScalculator.h"


// class ANT;
// class GScalculator;

class PBD_GSComponents
{
public:
  PBD_GSComponents(SimulationConfig* config);
  ~PBD_GSComponents();
  void CompoLogSetUp(Logger& logger);
  ANT* ant_;
  GScalculator* gscalculator_;

private:
};
