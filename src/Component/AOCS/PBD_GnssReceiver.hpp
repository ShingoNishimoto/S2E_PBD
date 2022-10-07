#pragma once

#include <Vector.hpp>
#include <GNSSReceiver.h>


using libra::Vector;

class PBD_GNSSReceiver : public GNSSReceiver
{
public:
  PBD_GNSSReceiver(const int prescaler, ClockGenerator* clock_gen, const int id, const std::string gnss_id, const int ch_max,
               const AntennaModel antenna_model, const Vector<3> ant_pos_b, const Quaternion q_b2c, const double half_width,
               const Vector<3> noise_std, const Vector<3> alignment_err_std, const Dynamics* dynamics, const GnssSatellites* gnss_satellites, const SimTime* simtime);

  void MainRoutine(int count);

  // Getter
  // Antenna Reference Point
  inline const Vector<3> GetAntennaPositionBody(void) const { return antenna_position_b_; }
  inline const Vector<3> GetAntennaPositionDesignECI(void) const { return arp_eci_.design_pos_; }
  inline const Vector<3> GetAntennaPositionTrueECI(void) const { return arp_eci_.true_pos_; }
  // Carrier Phase
  const Vector<3> GetPhaseReceivePositionTrueECI(void);
  inline const Vector<3> GetPhaseReceivePositionDesignECI(void) const { return phase_receive_position_eci_.design_pos_; }
  const Vector<3> GetPhaseReceivePositionDesignECI(const libra::Vector<3> sat_position) const;
  // inline const Vector<3> GetCPReceivePositionECEF(void) const { return phase_receive_position_ecef_; }
  // Pseudo Range
  const Vector<3> GetCodeReceivePositionTrueECI(void);
  inline const Vector<3> GetCodeReceivePositionDesignECI(void) const { return code_receive_position_eci_.design_pos_; }
  const Vector<3> GetCodeReceivePositionDesignECI(const libra::Vector<3> sat_position) const;
  // inline const Vector<3> GetPRReceivePositionECEF(void) const { return code_receive_position_ecef_; }
  // Alinment Error
  inline const Vector<3> GetAlignmentError(void) const { return alignment_err_b_; }

  std::string GetLogHeader() const;
  std::string GetLogValue() const;

protected:
  typedef struct
  {
    Vector<3> true_pos_{ 0.0 };
    Vector<3> design_pos_{ 0.0 };
  } PositionInfo;

  PositionInfo arp_eci_; // Antenna Reference Point position in ECI frame [m]
  PositionInfo phase_receive_position_eci_; // Antenna receive point in ECI frame [m]
  PositionInfo code_receive_position_eci_;  // Antenna receive point in ECI frame [m]
  Vector<3> alignment_err_b_{ 0.0 };          // Antenna alignment error

  libra::NormalRand nrs_antenna_b_x_, nrs_antenna_b_y_,
      nrs_antenna_b_z_;  // Random Error for each axis

  void UpdatePosition(void);
  void UpdateReceivePosition(Quaternion q_i2b);
};