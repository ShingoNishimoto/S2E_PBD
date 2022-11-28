#include "PBD_GroundStation.h"

PBD_GroundStation::PBD_GroundStation(SimulationConfig* config, int gs_id)
  :GroundStation(config, gs_id)
{
  Initialize(config);
}

PBD_GroundStation::~PBD_GroundStation()
{
  delete components_;
}

void PBD_GroundStation::Initialize(SimulationConfig* config)
{
  components_ = new PBD_GSComponents(config);
}

void PBD_GroundStation::LogSetup(Logger& logger)
{
  GroundStation::LogSetup(logger);
  components_->CompoLogSetUp(logger);
}

void PBD_GroundStation::Update(const Spacecraft& spacecraft, const GlobalEnvironment& global_env, const Antenna& sc_ant, const PBD_GroundStation& PBD_gs)
{
  GroundStation::Update(global_env.GetCelesInfo().GetEarthRotation(), spacecraft);
  components_->gscalculator_->Update(spacecraft, sc_ant, PBD_gs, *(components_->ant_));  // compo->ant_がnullの場合未定義動作になる気がするので対処が必要？
}
