#include "PBD_GNSSReceiver.hpp"
#include "../../Simulation/Spacecraft/PBD_Components.h"
#include <Library/math/GlobalRand.h>

#define FIXED_ALIGNMENT

PBD_GNSSReceiver::PBD_GNSSReceiver(const int prescaler, ClockGenerator* clock_gen, const int id, const std::string gnss_id, const int ch_max,
               const AntennaModel antenna_model, const Vector<3> ant_pos_b, const Quaternion q_b2c, const double half_width,
               const Vector<3> noise_std, const Vector<3> alignment_err_std,
               libra::Vector<3> pco, std::vector<double> pcv,
               const double azi_increment, const double ele_increment,
               const Dynamics* dynamics, const GnssSatellites* gnss_satellites, const SimTime* simtime):
               GNSSReceiver(prescaler, clock_gen, id, gnss_id, ch_max, antenna_model, ant_pos_b, q_b2c, half_width,
                            noise_std, dynamics, gnss_satellites, simtime),
               nrs_antenna_b_x_(0.0, alignment_err_std[0], g_rand.MakeSeed()),
               nrs_antenna_b_y_(0.0, alignment_err_std[1], g_rand.MakeSeed()),
               nrs_antenna_b_z_(0.0, alignment_err_std[2], g_rand.MakeSeed()),
               pcc_(PhaseCenterCorrection(pco, pcv, azi_increment, ele_increment))
{
#ifdef FIXED_ALIGNMENT
// stdではなくて，固定誤差として与えてあげる．
  alignment_err_b_[0] = alignment_err_std[0];
  alignment_err_b_[1] = alignment_err_std[1];
  alignment_err_b_[2] = alignment_err_std[2];
#else
  alignment_err_b_[0] = nrs_antenna_b_x_;
  alignment_err_b_[1] = nrs_antenna_b_y_;
  alignment_err_b_[2] = nrs_antenna_b_z_;
#endif // FIXED_ALIGNMENT

  int gnss_num = gnss_satellites_->GetNumOfSatellites();
  pre_observed_status_.assign(gnss_num, false);
  now_observed_status_.assign(gnss_num, false);

}

void PBD_GNSSReceiver::MainRoutine(int count)
{
  UNUSED(count);

  Vector<3> pos_true_eci_ = dynamics_->GetOrbit().GetSatPosition_i();
  Quaternion q_i2b = dynamics_->GetQuaternion_i2b();

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

void PBD_GNSSReceiver::CheckAntenna(const Vector<3> pos_true_eci_, Quaternion q_i2b) {
  if (antenna_model_ == SIMPLE)
    CheckAntennaSimple(pos_true_eci_, q_i2b);
  else if (antenna_model_ == CONE)
    CheckAntennaCone(pos_true_eci_, q_i2b);
}

void PBD_GNSSReceiver::CheckAntennaCone(const Vector<3> pos_true_eci_, Quaternion q_i2b) {
  // Cone model
  Vector<3> gnss_sat_pos_i, ant_pos_i, ant2gnss_i, ant2gnss_i_n, sat2ant_i;
  pre_observed_status_ = now_observed_status_; // copy
  vec_stocked_gnss_info_.clear();
  vec_gnssinfo_.clear();

  // antenna normal vector at inertial frame
  Vector<3> antenna_direction_c(0.0);
  antenna_direction_c[2] = 1.0;
  Vector<3> antenna_direction_b = q_b2c_.frame_conv_inv(antenna_direction_c);
  Vector<3> antenna_direction_i = q_i2b.frame_conv_inv(antenna_direction_b);

  sat2ant_i = q_i2b.frame_conv_inv(antenna_position_b_);
  ant_pos_i = pos_true_eci_ + sat2ant_i;

  // initialize
  gnss_sats_visible_num_ = 0;

  int gnss_num = gnss_satellites_->GetNumOfSatellites();
  now_observed_status_.assign(gnss_num, false); // numが変化するならpreとの整合性が取れない．．．

  for (int i = 0; i < gnss_num; i++) {
    // check if gnss ID is compatible with the receiver
    std::string id_tmp = gnss_satellites_->GetIDFromIndex(i);
    if (gnss_id_.find(id_tmp[0]) == std::string::npos) continue;

    // compute direction from sat to gnss in body-fixed frame
    gnss_sat_pos_i = gnss_satellites_->GetSatellitePositionEci(i);
    ant2gnss_i = gnss_sat_pos_i - ant_pos_i;
    double normalizer = 1 / norm(ant2gnss_i);
    ant2gnss_i_n = normalizer * ant2gnss_i;

    // check gnss sats are visible from antenna
    double Re = environment::earth_equatorial_radius_m;
    double inner1 = inner_product(ant_pos_i, gnss_sat_pos_i);
    int is_visible_ant2gnss = 0;
    if (inner1 > 0)
      is_visible_ant2gnss = 1;
    else {
      Vector<3> tmp = ant_pos_i + inner_product(-ant_pos_i, ant2gnss_i_n) * ant2gnss_i;
      if (norm(tmp) < Re)
        // There is earth between antenna and gnss
        is_visible_ant2gnss = 0;
      else
        // There is not earth between antenna and gnss
        is_visible_ant2gnss = 1;
    }

    double inner2 = inner_product(antenna_direction_i, ant2gnss_i_n);
    if (inner2 > cos(half_width_ * libra::deg_to_rad) && is_visible_ant2gnss) {
      // is visible
      if (pre_observed_status_.at(i)) // 見えていたものから優先的に
      {
        gnss_sats_visible_num_++;
        SetGnssInfo(ant2gnss_i, q_i2b, id_tmp);
        now_observed_status_.at(i) = true;
      }
      else
      {
        SetStockedGnssInfo(ant2gnss_i, q_i2b, id_tmp);
      }
    }
    if (gnss_sats_visible_num_ >= ch_max_) break; // 一旦ch以上は観測しないようにbreakする．もう少しちゃんとした選択アルゴリズムを実装する．
  }

  if (gnss_sats_visible_num_ < ch_max_)
  {
    // stockされているものを上から追加する．
    for (int i = 0; i < vec_stocked_gnss_info_.size(); i++)
    {
      const GnssInfo gnss_info = vec_stocked_gnss_info_.at(i);
      const int id = gnss_satellites_->GetIndexFromID(gnss_info.ID);
      vec_gnssinfo_.push_back(gnss_info);
      now_observed_status_.at(id) = true;
      gnss_sats_visible_num_++;
      if (gnss_sats_visible_num_ >= ch_max_) break;
    }
  }

  if (gnss_sats_visible_num_ > 0)
    is_gnss_sats_visible_ = 1;
  else
    is_gnss_sats_visible_ = 0;
}


void PBD_GNSSReceiver::SetStockedGnssInfo(Vector<3> ant2gnss_i, Quaternion q_i2b, std::string gnss_id) {
  Vector<3> ant2gnss_b, ant2gnss_c;

  ant2gnss_b = q_i2b.frame_conv(ant2gnss_i);
  ant2gnss_c = q_b2c_.frame_conv(ant2gnss_b);

  double dist = norm(ant2gnss_c);
  double lon = AcTan(ant2gnss_c[1], ant2gnss_c[0]);
  double lat = AcTan(ant2gnss_c[2], sqrt(pow(ant2gnss_c[0], 2.0) + pow(ant2gnss_c[1], 2.0)));

  GnssInfo gnss_info_new = {gnss_id, lat, lon, dist};
  vec_stocked_gnss_info_.push_back(gnss_info_new);
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
  Vector<3> antenna_direction_c = (1 / libra::norm(antenna_position_b_)) * antenna_position_b_;
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
