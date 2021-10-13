#include "PBD_GSComponents.h"
#include "Initialize.h"


PBD_GSComponents::PBD_GSComponents(SimulationConfig* config)
{
  IniAccess iniAccess = IniAccess(config->ini_base_fname_);
	
  string ant_ini_path = iniAccess.ReadString("COMPONENTS_FILE", "ant_gs_file");
  config->main_logger_->CopyFileToLogDir(ant_ini_path);
  ant_ = new ANT(InitANT(1, ant_ini_path));		//TODO : deal with multiple ANT
  string gscalculator_ini_path = iniAccess.ReadString("COMPONENTS_FILE", "gscalculator_file");
  config->main_logger_->CopyFileToLogDir(gscalculator_ini_path);
  gscalculator_ = new GScalculator(InitGScalculator(gscalculator_ini_path));  // GScalcはGSごとに固有のものなのでidは不要か
}

PBD_GSComponents::~PBD_GSComponents()
{
  delete ant_;
  delete gscalculator_;
}

void PBD_GSComponents::CompoLogSetUp(Logger& logger)
{
  // logger.AddLoggable(ant_);
  logger.AddLoggable(gscalculator_);
}
