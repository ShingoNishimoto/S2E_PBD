#include "PBD_GNSSReceiver.hpp"
#include "../../Simulation/Spacecraft/PBD_Components.h"
#include <Library/math/GlobalRand.h>

PBD_GNSSReceiver::PBD_GNSSReceiver(const int prescaler, ClockGenerator* clock_gen, const int id, const std::string gnss_id, const int ch_max,
               const AntennaModel antenna_model, const Vector<3> ant_pos_b, const Quaternion q_b2c, const double half_width,
               const Vector<3> noise_std, const Vector<3> alignment_err_std, const Dynamics* dynamics, const GnssSatellites* gnss_satellites, const SimTime* simtime):
               GNSSReceiver(prescaler, clock_gen, id, gnss_id, ch_max, antenna_model, ant_pos_b, q_b2c, half_width,
                            noise_std, dynamics, gnss_satellites, simtime),
               nrs_antenna_b_x_(0.0, alignment_err_std[0], g_rand.MakeSeed()),
               nrs_antenna_b_y_(0.0, alignment_err_std[1], g_rand.MakeSeed()),
               nrs_antenna_b_z_(0.0, alignment_err_std[2], g_rand.MakeSeed())
{
  alignment_err_b_[0] = nrs_antenna_b_x_;
  alignment_err_b_[1] = nrs_antenna_b_y_;
  alignment_err_b_[2] = nrs_antenna_b_z_;
}

void PBD_GNSSReceiver::MainRoutine(int count)
{
  UNUSED(count);

  Vector<3> pos_true_eci_ = dynamics_->GetOrbit().GetSatPosition_i();
  Quaternion q_i2b = dynamics_->GetQuaternion_i2b();

  // 特定のGNSS衛星に対してみているのではなくて，反地球方向を見ているかの簡単な確認．
  CheckAntenna(pos_true_eci_, q_i2b);

  if (is_gnss_sats_visible_ == 1) {  // Antenna of GNSS-R can detect GNSS signal
    UpdatePosition();
    // AddNoise(pos_true_eci_, position_ecef_); noiseをこっちでは追加しない．
    UpdateReceivePosition(q_i2b);

    utc_ = simtime_->GetCurrentUTC();
    ConvertJulianDayToGPSTime(simtime_->GetCurrentJd());
  } else {
    // position information will not be updated in this case
    // only time information will be updated in this case (according to the
    // receiver's internal clock)
    utc_ = simtime_->GetCurrentUTC();
    ConvertJulianDayToGPSTime(simtime_->GetCurrentJd());
  }
}

void PBD_GNSSReceiver::UpdatePosition(void)
{
  position_ecef_ = dynamics_->GetOrbit().GetSatPosition_ecef();
  position_eci_ = dynamics_->GetOrbit().GetSatPosition_i();
  velocity_ecef_ = dynamics_->GetOrbit().GetSatVelocity_ecef();
  velocity_eci_ = dynamics_->GetOrbit().GetSatVelocity_i();
  position_llh_ = dynamics_->GetOrbit().GetLatLonAlt();
}

void PBD_GNSSReceiver::UpdateReceivePosition(Quaternion q_i2b)
{
  Vector<3> gnss_sat_pos_i, sat2ant_design_i, sat2ant_true_i;

  // antenna normal vector at inertial frame
  Vector<3> antenna_direction_c(0.0);
  antenna_direction_c[2] = 1.0; // PZに搭載
  Vector<3> antenna_direction_b = q_b2c_.frame_conv_inv(antenna_direction_c);
  Vector<3> antenna_direction_i = q_i2b.frame_conv_inv(antenna_direction_b);

  sat2ant_design_i = q_i2b.frame_conv_inv(antenna_position_b_);
  sat2ant_true_i = sat2ant_design_i + q_i2b.frame_conv_inv(alignment_err_b_);
  arp_eci_.design_pos_ = position_eci_ + sat2ant_design_i;
  arp_eci_.true_pos_ = position_eci_ + sat2ant_true_i;

  // code position
  code_receive_position_eci_ = arp_eci_; // ひとまず一致させる．
  // phase position
  phase_receive_position_eci_ = arp_eci_; // TODO: PCOを考える．

  // とりあえずecefは置いておく．
  // Matrix<3, 3> trans_eci2ecef_ = local_env.GetCelesInfo().GetGlobalInfo().GetEarthRotation().GetDCMJ2000toXCXF();
}

// ここはGNSS受信機が出すタイミングに合わせるのが普通なので，このタイミングで更新をするのはやらないが，一旦この方針，TODO: 後々修正
const Vector<3> PBD_GNSSReceiver::GetPhaseReceivePositionTrueECI(void)
{
  UpdatePosition();
  Quaternion q_i2b = dynamics_->GetQuaternion_i2b();
  UpdateReceivePosition(q_i2b);

  return phase_receive_position_eci_.true_pos_;
}

const Vector<3> PBD_GNSSReceiver::GetCodeReceivePositionTrueECI(void)
{
  UpdatePosition();
  Quaternion q_i2b = dynamics_->GetQuaternion_i2b();
  UpdateReceivePosition(q_i2b);

  return code_receive_position_eci_.true_pos_;
}

const Vector<3> PBD_GNSSReceiver::GetPhaseReceivePositionDesignECI(const libra::Vector<3> sat_position) const
{
  Vector<3> sat2ant_design_i, receive_position_i;
  Quaternion q_i2b = dynamics_->GetQuaternion_i2b();

  sat2ant_design_i = q_i2b.frame_conv_inv(antenna_position_b_);
  receive_position_i = sat_position + sat2ant_design_i; // TODO: add pco, pcv

  return receive_position_i;
}

const Vector<3> PBD_GNSSReceiver::GetCodeReceivePositionDesignECI(const libra::Vector<3> sat_position) const
{
  Vector<3> sat2ant_design_i, receive_position_i;
  Quaternion q_i2b = dynamics_->GetQuaternion_i2b();

  sat2ant_design_i = q_i2b.frame_conv_inv(antenna_position_b_);
  receive_position_i = sat_position + sat2ant_design_i; // TODO: add pco, pcv

  return receive_position_i;
}


std::string PBD_GNSSReceiver::GetLogHeader() const  // For logs
{
  std::string str_tmp = "";
  /*
  str_tmp += WriteScalar("gnss_year");
  str_tmp += WriteScalar("gnss_month");
  str_tmp += WriteScalar("gnss_day");
  str_tmp += WriteScalar("gnss_hour");
  str_tmp += WriteScalar("gnss_min");
  str_tmp += WriteScalar("gnss_sec");
  */
  str_tmp += WriteVector("gnss_position", "eci", "m", 3);
  str_tmp += WriteVector("gnss_arp_true", "eci", "m", 3);
  str_tmp += WriteVector("gnss_arp_design", "eci", "m", 3);
  str_tmp += WriteVector("gnss_velocity", "ecef", "m/s", 3);
  str_tmp += WriteScalar("gnss_lat", "rad");
  str_tmp += WriteScalar("gnss_lon", "rad");
  str_tmp += WriteScalar("gnss_alt", "m");
  str_tmp += WriteScalar("gnss_vis_flag");
  str_tmp += WriteScalar("gnss_vis_num");

  return str_tmp;
}

std::string PBD_GNSSReceiver::GetLogValue() const  // For logs
{
  std::string str_tmp = "";
  /*
  str_tmp += WriteScalar(utc_.year);
  str_tmp += WriteScalar(utc_.month);
  str_tmp += WriteScalar(utc_.day);
  str_tmp += WriteScalar(utc_.hour);
  str_tmp += WriteScalar(utc_.min);
  str_tmp += WriteScalar(utc_.sec);
  */
  str_tmp += WriteVector(position_eci_, 10);
  str_tmp += WriteVector(arp_eci_.true_pos_, 10);
  str_tmp += WriteVector(arp_eci_.design_pos_, 10);
  str_tmp += WriteVector(velocity_ecef_, 10);
  str_tmp += WriteScalar(position_llh_[0], 10);
  str_tmp += WriteScalar(position_llh_[1], 10);
  str_tmp += WriteScalar(position_llh_[2], 10);
  str_tmp += WriteScalar(is_gnss_sats_visible_);
  str_tmp += WriteScalar(gnss_sats_visible_num_);

  return str_tmp;
}
