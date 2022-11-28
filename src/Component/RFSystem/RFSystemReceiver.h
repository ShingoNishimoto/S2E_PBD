#pragma once
#include "../../Simulation/InterSatComm/RFSystem/RFSystemBeam.h"
#include "../../Simulation/InterSatComm/PBD_InterSatComm.h"
#include "Antenna.hpp"
#include "Vector.hpp"
#include "ComponentBase.h"
#include "ILoggable.h"
#include "Dynamics.h"
#include "Quaternion.hpp"
#include <utility>

class RFSystemReceiver : public ComponentBase, public ILoggable
{
public:
  RFSystemReceiver(
    const int prescaler,
    ClockGenerator* clock_gen,
    const Vector<3>& compo_position_b,
    const Quaternion& q_b2c,
    Antenna* ant,
    PBD_InterSatComm* pbd_inter_sat_comm,
    const Dynamics* dynamics,
    const double update_interval_sec);
  ~RFSystemReceiver();

  void MainRoutine(int count) override;
  virtual std::string GetLogHeader() const;
  virtual std::string GetLogValue() const;

  //Getter
  inline double GetReceivedPower() const { return received_power_watt_; }

private:
  //Component Parameters
  Vector<3> compo_position_b_; //Position in relation to the center of gravity of the satellite[m]
  Vector<3> line_of_sight_c_;
  Quaternion q_b2c_;
  double update_interval_sec_;

  //State Values
  double received_power_watt_ = 0.0; //[W]

  enum class RFLaserLinkState{
    Lost,
    WaitingForLink,
    Linked
  };
  RFLaserLinkState state_ = RFLaserLinkState::Lost;
  double link_wait_counter_sec_ = 0.0;

  //Menber Functions
  double CalcReceivedPower();
  void UpdateLinkState();

  //External Data
  const Dynamics* dynamics_;
  RFSystemBeam* rf_system_beam_;
  Antenna* ant_;
};
