#include "../Initialize_PBD.h"

RFSystemTransmitter InitRFSystemTransmitter(ClockGenerator* clock_gen, const std::string ini_path, PBD_InterSatComm* pbd_inter_sat_comm, const Dynamics* dynamics)
{
  auto conf = IniAccess(ini_path);
  char* section = "RFSystemTransmitter";

  int prescaler = conf.ReadInt(section, "prescaler");

  Vector<3> compo_position_b;
  conf.ReadVector(section, "compo_position_b", compo_position_b);

  Quaternion q_b2c;
  conf.ReadQuaternion(section, "q_b2c", q_b2c);

  q_b2c = q_b2c.normalize();

  std::string ant_ini_path = conf.ReadString("COMPONENT_FILE", "rf_system_transmitter_file");
  Antenna ant_ = Antenna(InitAntenna(2, ant_ini_path));

  RFSystemTransmitter rf_sys_t(
    prescaler,
    clock_gen,
    compo_position_b,
    q_b2c,
    &ant_,
    pbd_inter_sat_comm,
    dynamics
  );

  return rf_sys_t;
}

