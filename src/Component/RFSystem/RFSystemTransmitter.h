#pragma once
#include "../../Simulation/InterSatComm/RFSystem/RFSystemBeam.h"
#include "Antenna.hpp"
#include "../../Simulation/InterSatComm/PBD_InterSatComm.h"
#include "ComponentBase.h"
#include "ILoggable.h"
#include "Dynamics.h"
#include "Quaternion.hpp"
#include "Vector.hpp"

class RFSystemTransmitter : public ComponentBase, public ILoggable
{
public:
  RFSystemTransmitter(
    int prescaler,
    ClockGenerator* clock_gen,
    const Vector<3>& compo_position_b,
    const Quaternion& q_b2c,
    Antenna* ant,
    PBD_InterSatComm* pbd_inter_sat_comm,
    const Dynamics* dynamics
  );

  ~RFSystemTransmitter();

  void MainRoutine(int count) override;
  virtual std::string GetLogHeader() const;
  virtual std::string GetLogValue() const;

  //setter
  // inline void InputAcousticSignal(Vector<2> acoustic_signal_freq_hz) { acoustic_signal_freq_hz_ = acoustic_signal_freq_hz; }

private:

  Vector<3> compo_position_b_{ 0.0 };//Position in relation to the center of gravity of the satellite[m]
  Vector<3> line_of_sight_c_;//[ ]
  Quaternion q_b2c_;

  const Dynamics* dynamics_;
  RFSystemBeam* rf_system_beam_;
  Antenna* ant_;

  //private functions
  //Vector<3> CalcPointingVector_c();
  void ExpandBeam();
  void AddNoise();
};
