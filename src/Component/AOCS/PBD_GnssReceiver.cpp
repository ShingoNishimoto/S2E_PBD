#include "PBD_GNSSReceiver.hpp"
#include "../../Simulation/Spacecraft/PBD_Components.h"
#include <Library/math/GlobalRand.h>

PBD_GnssReceiver::PBD_GnssReceiver(const int prescaler, ClockGenerator* clock_gen, const int id, const std::string gnss_id, const int ch_max,
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

void PBD_GnssReceiver::MainRoutine(int count)
{
  UNUSED(count);

  Vector<3> pos_true_eci_ = dynamics_->GetOrbit().GetSatPosition_i();
  Quaternion q_i2b = dynamics_->GetQuaternion_i2b();

  // 特定のGNSS衛星に対してみているのではなくて，反地球方向を見ているかの簡単な確認．
  CheckAntenna(pos_true_eci_, q_i2b); // simple_modelにしておく．

  if (is_gnss_sats_visible_ == 1) {  // Antenna of GNSS-R can detect GNSS signal
    position_ecef_ = dynamics_->GetOrbit().GetSatPosition_ecef();
    position_eci_ = dynamics_->GetOrbit().GetSatPosition_i();
    velocity_ecef_ = dynamics_->GetOrbit().GetSatVelocity_ecef();
    velocity_eci_ = dynamics_->GetOrbit().GetSatVelocity_i();
    position_llh_ = dynamics_->GetOrbit().GetLatLonAlt();
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

// ここはCheckAntennaの部分をオーバーライドしてもいいかも
void PBD_GnssReceiver::UpdateReceivePosition(Quaternion q_i2b)
{
  Vector<3> gnss_sat_pos_i, sat2ant_i;

  // antenna normal vector at inertial frame
  Vector<3> antenna_direction_c(0.0);
  antenna_direction_c[2] = 1.0;
  Vector<3> antenna_direction_b = q_b2c_.frame_conv_inv(antenna_direction_c);
  Vector<3> antenna_direction_i = q_i2b.frame_conv_inv(antenna_direction_b);

  sat2ant_i = q_i2b.frame_conv_inv(antenna_position_b_ + alignment_err_b_);
  antenna_position_eci_ = position_eci_ + sat2ant_i;

  // code position
  code_receive_point_eci_ = antenna_position_eci_; // ひとまず一致させる．
  // phase position
  phase_receive_point_eci_ = antenna_position_eci_; // TODO: PCOを考える．

  // とりあえずecefは置いておく．
  // Matrix<3, 3> trans_eci2ecef_ = local_env.GetCelesInfo().GetGlobalInfo().GetEarthRotation().GetDCMJ2000toXCXF();
}
