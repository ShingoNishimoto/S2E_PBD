#include "../Initialize_PBD.h"
#include "../../../Simulation/InterSatComm/RFSystem/RFSystemBeam.h"

RFSystemBeam InitRFSystemBeam(const string ini_path)
{
  auto conf = IniAccess(ini_path);
  char* section = "RFSystemBeam";

  double wavelength_m = conf.ReadDouble(section, "wavelength");
  double r_beam_waist_m = conf.ReadDouble(section, "r_beam_waist");
  double total_power_watt = conf.ReadDouble(section, "total_power");


  string name_tag = conf.ReadString(section, "name_tag");

  RFSystemBeam rf_system_beam;

  return rf_system_beam;
}
