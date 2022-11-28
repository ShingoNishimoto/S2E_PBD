#pragma once
#include <Component/CommGS/InitAntenna.hpp>
#include <Component/CommGS/InitGsCalculator.hpp>

class PBD_GSComponents
{
public:
  PBD_GSComponents(SimulationConfig* config);
  ~PBD_GSComponents();
  void CompoLogSetUp(Logger& logger);
  Antenna* ant_;
  GScalculator* gscalculator_;

private:
};
