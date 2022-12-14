#pragma once

#include <Vector.hpp>
#include <GNSSReceiver.h>
#include "../../Simulation/PBD/PhaseCenterCorrection.hpp"

using libra::Vector;

class PBD_GNSSReceiver : public GNSSReceiver
{
public:
  PBD_GNSSReceiver(const int prescaler, ClockGenerator* clock_gen, const int id,
    const std::string gnss_id, const int ch_max, const AntennaModel antenna_model,
    const Vector<3> ant_pos_b, const Quaternion q_b2c, const double half_width,
    const Vector<3> noise_std, const Vector<3> alignment_err_std,
    libra::Vector<3> pco, std::vector<double> pcv,
    const double azi_increment, const double ele_increment,
    const Dynamics* dynamics, const GnssSatellites* gnss_satellites,
    const SimTime* simtime);

  void MainRoutine(int count);

  // Getter
  // Antenna Reference Point
  inline const Vector<3> GetAntennaPositionBody(void) const { return antenna_position_b_; }
  inline const Vector<3> GetAntennaPositionDesignECI(void) const { return arp_eci_.design_pos_; }
  inline const Vector<3> GetAntennaPositionTrueECI(void) const { return arp_eci_.true_pos_; }
  // Carrier Phase
  const Vector<3> GetPhaseReceivePositionTrueECI(void);
  // inline const Vector<3> GetPhaseReceivePositionDesignECI(void) const {return phase_receive_position_eci_.design_pos_;}
  const Vector<3> GetPhaseReceivePositionDesignECI(const libra::Vector<3> sat_position) const;
  // inline const Vector<3> GetPhaseReceivePositionECEF(void) const { return phase_receive_position_ecef_; }
  // Pseudo Range
  const Vector<3> GetCodeReceivePositionTrueECI(void);
  // inline const Vector<3> GetCodeReceivePositionDesignECI(void) const { return code_receive_position_eci_.design_pos_; }
  const Vector<3> GetCodeReceivePositionDesignECI(const libra::Vector<3> sat_position) const;
  // inline const Vector<3> GetCodeReceivePositionECEF(void) const { return code_receive_position_ecef_; }
  libra::Vector<3> TransCompoToEci(const libra::Vector<3>& target_vec_c);

  // Alignment Error
  inline const Vector<3> GetAlignmentError(void) const { return alignment_err_b_; }
  inline const double GetPCC_m(const double azimuth_deg, const double elevation_deg) { return pcc_.GetPCC_m(azimuth_deg, elevation_deg); }
  inline const libra::Vector<3> GetPCO_mm(void) const { return pcc_.GetPCO_mm(); }
  inline const vector<GnssInfo> GetGnssInfoVec(void) const { return vec_gnssinfo_; }
  inline const double GetClockBias(void) const { return clock_bias_; }

  struct GnssReceiverObservations
  {
    double l1_pseudo_range;
    double l2_pseudo_range;
    std::pair<double, double> l1_carrier_phase; //[位相, N(整数, dtype = double)]
    std::pair<double, double> l2_carrier_phase; //[位相, N(整数, dtype = double)]
  };

  GnssReceiverObservations GetRawObservations(const int ch);

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
  PhaseCenterCorrection pcc_;

  libra::NormalRand nrs_antenna_b_x_, nrs_antenna_b_y_,
      nrs_antenna_b_z_;  // Random Error for each axis
  std::vector<bool> pre_observed_status_{};
  std::vector<bool> now_observed_status_{};
  std::vector<GnssInfo> vec_stocked_gnss_info_{};
  double clock_bias_; // receiver clock biasの真値[m]
  // std::random_device seed_gen;
  std::mt19937 mt_;

  void UpdatePosition(void);
  void UpdateReceivePosition(Quaternion q_i2b);
  void CheckAntenna(const Vector<3> pos_true_eci_, Quaternion q_i2b);
  void CheckAntennaCone(const Vector<3> pos_true_eci_, Quaternion q_i2b);
  void SetStockedGnssInfo(Vector<3> ant2gnss_i, Quaternion q_i2b, std::string gnss_id);
};
