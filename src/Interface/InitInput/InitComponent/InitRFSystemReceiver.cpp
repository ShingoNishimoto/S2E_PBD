#include "../Initialize_PBD.h"

RFSystemReceiver InitRFSystemReceiver(ClockGenerator* clock_gen, const std::string ini_path, PBD_InterSatComm* pbd_inter_sat_comm, const Dynamics* dynamics, const SimTime* sim_time)
{
  auto conf = IniAccess(ini_path);
  char* section = "RFSystemReceiver";

  int prescaler = conf.ReadInt(section, "prescaler");

  Vector<3> compo_position_b;
  conf.ReadVector(section, "compo_position_b", compo_position_b);

  Quaternion q_b2c;
  conf.ReadQuaternion(section, "q_b2c", q_b2c);

  double update_interval_sec = prescaler * sim_time->GetCompoStepSec();

  std::string ant_ini_path = conf.ReadString("COMPONENT_FILE", "rf_system_receiver_file");
  Antenna ant_ = Antenna(InitAntenna(1, ant_ini_path));

  RFSystemReceiver rf_sys_r(
    prescaler,
    clock_gen,
    compo_position_b,
    q_b2c,
    &ant_,
    pbd_inter_sat_comm,
    dynamics,
    update_interval_sec
  );

  return rf_sys_r;
}

