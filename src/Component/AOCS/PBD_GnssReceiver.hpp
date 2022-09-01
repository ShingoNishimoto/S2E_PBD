#pragma once

#include <Vector.hpp>
#include <GNSSReceiver.h>


using libra::Vector;

class PBD_GNSSReceiver : public GNSSReceiver
{
public:
  PBD_GNSSReceiver(const int prescaler, ClockGenerator* clock_gen, const int id, const std::string gnss_id, const int ch_max,
               const AntennaModel antenna_model, const Vector<3> ant_pos_b, const Quaternion q_b2c, const double half_width,
               const Vector<3> noise_std, const Vector<3> alignment_err_std, const Dynamics* dynamics, const GNSSSatellites* gnss_satellites, const SimTime* simtime);

  void MainRoutine(int count);

  // Getter
  inline const Vector<3> GetCPReceivePositionECI(void) { return phase_receive_point_eci_; }
  inline const Vector<3> GetCPReceivePositionECEF(void) { return phase_receive_point_ecef_; }
  inline const Vector<3> GetPRReceivePositionECI(void) { return code_receive_point_eci_; }
  inline const Vector<3> GetPRReceivePositionECEF(void) { return code_receive_point_ecef_; }

protected:
  Vector<3> alignment_err_b_; // Antenna alignment error
  Vector<3> antenna_position_eci_; // Antenna Reference Point position in ECI frame
  Vector<3> phase_receive_point_eci_; // Antenna receive point in ECI frame<- PCVを考えない簡易的なもの．
  Vector<3> phase_receive_point_ecef_; // Antenna receive point in ECEF frame
  Vector<3> code_receive_point_eci_;  // Antenna receive point in ECI frame <- これは普通にARPに等しいと考えていい？
  Vector<3> code_receive_point_ecef_;  // Antenna receive point in ECEF frame
  libra::NormalRand nrs_antenna_b_x_, nrs_antenna_b_y_,
      nrs_antenna_b_z_;  // Random Error for each axis
  
  void UpdateReceivePosition(Quaternion q_i2b);
};
