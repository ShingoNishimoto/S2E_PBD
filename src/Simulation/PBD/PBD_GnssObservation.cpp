#include "PBD_const.h"
#include "PBD_GnssObservation.h"

bool GnssObservedValues::check_normal()
{
    unsigned int n = observable_gnss_sat_id.size();
    if(n != gnss_satellites_position.size()) return false;
    if(n != gnss_clock.size()) return false;
    if(n != L1_pseudo_range.size()) return false;
    if(n != L2_pseudo_range.size()) return false;
    if(n != L1_carrier_phase.size()) return false;
    if(n != L2_carrier_phase.size()) return false;
    if(n != ionfree_pseudo_range.size()) return false;
    //if(n != ionfree_carrier_phase.size()) return false;

    return true;
}

PBD_GnssObservation::PBD_GnssObservation(const Orbit& orbit, const GnssSatellites& gnss_satellites) : orbit_(orbit), gnss_satellites_(gnss_satellites)
{
  num_of_gnss_satellites_ = gnss_satellites_.GetNumOfSatellites();
  // 初期化
  info_.pre_observed_status.assign(num_of_gnss_satellites_, false);
  info_.now_observed_status.assign(num_of_gnss_satellites_, false);
  std::normal_distribution<> receiver_clock_dist(0.0, clock_sigma);
  receiver_clock_bias_ = receiver_clock_dist(mt);
}

PBD_GnssObservation::~PBD_GnssObservation() {}

void PBD_GnssObservation::Update(void)
{
  // receiver clock
  std::normal_distribution<> receiver_clock_dist(0.0, clock_sigma);
  receiver_clock_bias_ = receiver_clock_dist(mt);
  // ここで毎回呼ぶと軌道情報が更新されてなくて死ぬ．
  UpdateGnssObservation();
}

void PBD_GnssObservation::UpdateGnssObservation()
{
  //推定値の計算
  num_of_gnss_satellites_ = gnss_satellites_.GetNumOfSatellites(); // 更新 <- ?
  libra::Vector<3> sat_position_i = orbit_.GetSatPosition_i(); // ECI

  info_.now_observed_status.assign(num_of_gnss_satellites_, false);
  info_.now_observed_gnss_sat_id.clear(); //クリア
  ClearPreValues(true_values_);
  ClearPreValues(observed_values_);

  for (int gnss_sat_id = 0; gnss_sat_id < num_of_gnss_satellites_; ++gnss_sat_id) {
    //if(gnss_sat_id == 7 || gnss_sat_id == 23 || gnss_sat_id == 31) continue; 　←この衛星たちの軌道情報が悪いからこうしていたのか？
    if (!gnss_satellites_.GetWhetherValid(gnss_sat_id)) continue;
    libra::Vector<3> gnss_position = gnss_satellites_.Get_true_info().GetSatellitePositionEci(gnss_sat_id);
    bool see_flag = CheckCanSeeSatellite(sat_position_i, gnss_position);

    if (!see_flag)
    {
      continue;
    }
    info_.now_observed_status.at(gnss_sat_id) = true;
    info_.now_observed_gnss_sat_id.push_back(gnss_sat_id);

    double gnss_clock = gnss_satellites_.Get_true_info().GetSatelliteClock(gnss_sat_id); // これはclock bias

    // 以下で衛星位置を与えたら重心間の距離として観測量が出てくる．この時点でアンテナ取付位置の慣性座標を渡す必要がある．
    double l1_pseudo_range = gnss_satellites_.GetPseudoRangeECI(gnss_sat_id, sat_position_i, receiver_clock_bias_, L1_frequency);
    double l2_pseudo_range = gnss_satellites_.GetPseudoRangeECI(gnss_sat_id, sat_position_i, receiver_clock_bias_, L2_frequency);
    // この中に整数不定性を入れてないのがダメなのでは？
    auto l1_carrier_phase = gnss_satellites_.GetCarrierPhaseECI(gnss_sat_id, sat_position_i, receiver_clock_bias_, L1_frequency);
    auto l2_carrier_phase = gnss_satellites_.GetCarrierPhaseECI(gnss_sat_id, sat_position_i, receiver_clock_bias_, L2_frequency);

    double ionfree_range = (pow(L1_frequency / L2_frequency, 2.0) * l1_pseudo_range - l2_pseudo_range) / (pow(L1_frequency / L2_frequency, 2.0) - 1);
    double ionfree_phase = (pow(L1_frequency / L2_frequency, 2.0) * L1_lambda * (l1_carrier_phase.first + l1_carrier_phase.second) - L2_lambda * (l2_carrier_phase.first + l2_carrier_phase.second)) / (pow(L1_frequency / L2_frequency, 2.0) - 1);

    true_values_.observable_gnss_sat_id.push_back(gnss_sat_id);
    true_values_.gnss_satellites_position.push_back(gnss_position);
    true_values_.gnss_clock.push_back(gnss_clock);
    true_values_.L1_pseudo_range.push_back(l1_pseudo_range);
    true_values_.L2_pseudo_range.push_back(l2_pseudo_range);
    true_values_.L1_carrier_phase.push_back(l1_carrier_phase);
    true_values_.L2_carrier_phase.push_back(l2_carrier_phase);

    true_values_.ionfree_pseudo_range.push_back(ionfree_range);
    true_values_.ionfree_carrier_phase.push_back(ionfree_phase);

    // 観測情報の方には観測誤差を混ぜる
    std::normal_distribution<> pseudo_range_noise(0.0, pseudo_sigma);
    std::normal_distribution<> carrier_phase_noise(0.0, carrier_sigma);

    //estimateに使う方の情報
    gnss_position = gnss_satellites_.GetSatellitePositionEci(gnss_sat_id);
    gnss_clock = gnss_satellites_.GetSatelliteClock(gnss_sat_id);

    // add measurement error
    l1_pseudo_range += pseudo_range_noise(mt);
    l2_pseudo_range += pseudo_range_noise(mt);
    l1_carrier_phase.first += carrier_phase_noise(mt) / L1_lambda;
    l2_carrier_phase.first += carrier_phase_noise(mt) / L2_lambda;

    ionfree_range = (pow(L1_frequency / L2_frequency, 2.0) * l1_pseudo_range - l2_pseudo_range) / (pow(L1_frequency / L2_frequency, 2.0) - 1);
    ionfree_phase = L2_lambda * (L1_frequency / L2_frequency * (l1_carrier_phase.first + l1_carrier_phase.second) - (l2_carrier_phase.first + l2_carrier_phase.second)) / (pow(L1_frequency / L2_frequency, 2.0) - 1);

    observed_values_.observable_gnss_sat_id.push_back(gnss_sat_id);
    observed_values_.gnss_satellites_position.push_back(gnss_position);
    observed_values_.gnss_clock.push_back(gnss_clock);
    observed_values_.L1_pseudo_range.push_back(l1_pseudo_range);
    observed_values_.L2_pseudo_range.push_back(l2_pseudo_range);
    observed_values_.L1_carrier_phase.push_back(l1_carrier_phase);
    observed_values_.L2_carrier_phase.push_back(l2_carrier_phase);

    observed_values_.ionfree_pseudo_range.push_back(ionfree_range);
    observed_values_.ionfree_carrier_phase.push_back(ionfree_phase);
  }
}


bool PBD_GnssObservation::CheckCanSeeSatellite(const libra::Vector<3> satellite_position, const libra::Vector<3> gnss_position) const
{
  // ここには姿勢の要素も入れなければいけない．
  double angle_rad = angle(satellite_position, gnss_position - satellite_position);
  if (angle_rad < M_PI / 2.0 - mask_angle) return true;
  else return false;
}

// ionfreeの計算をしている．
void PBD_GnssObservation::CalcIonfreeObservation()
{
  int observed_gnss_index = 0;
  for (int i = 0; i < num_of_gnss_satellites_; ++i)
  {
    if (info_.pre_observed_status.at(i) == true && info_.now_observed_status.at(i) == false)
    {
      l1_bias_.at(i) = 0.0;
      l2_bias_.at(i) = 0.0;
    }
    else if (info_.pre_observed_status.at(i) == false && info_.now_observed_status.at(i) == true)
    {
      // (first + second)*lambda から真の距離引いてそこからN求める．ここから！！！！！！！！！！！！！！！！！！！！！ そのまま真の距離引いたら0になるからここでの真の距離は時刻を使う．時刻の精度以下に埋もれる部分が整数不定性として出てくる？伝搬時間も必要やん．
      l1_bias_.at(i) = true_values_.L1_carrier_phase.at(observed_gnss_index).second; // これじゃダメ．あまり分のNを求めないと．あんま関係ない気がするので後で対応する．
      // これがどのchに対応しているかはわかっている．
      l2_bias_.at(i) = true_values_.L2_carrier_phase.at(observed_gnss_index).second;
    }
    if (info_.now_observed_status.at(i)) ++observed_gnss_index;
  }

  for (int i = 0; i < true_values_.observable_gnss_sat_id.size(); ++i)
  {
    int gnss_sat_id = true_values_.observable_gnss_sat_id.at(i);

    auto L1_observed = true_values_.L1_carrier_phase.at(i);
    L1_observed.first += L1_observed.second - l1_bias_.at(gnss_sat_id);

    auto L2_observed = true_values_.L2_carrier_phase.at(i);
    L2_observed.first += L2_observed.second - l2_bias_.at(gnss_sat_id);

    double ionfree_phase = (pow(L1_frequency / L2_frequency, 2.0) * L1_lambda * L1_observed.first - L2_lambda * L2_observed.first) / (pow(L1_frequency / L2_frequency, 2.0) - 1);
    true_values_.ionfree_carrier_phase.push_back(ionfree_phase);
  }

  for (int i = 0; i < observed_values_.observable_gnss_sat_id.size(); ++i)
  {
    int gnss_sat_id = observed_values_.observable_gnss_sat_id.at(i);

    auto L1_observed = observed_values_.L1_carrier_phase.at(i);
    L1_observed.first += L1_observed.second - l1_bias_.at(gnss_sat_id);

    auto L2_observed = observed_values_.L2_carrier_phase.at(i);
    L2_observed.first += L2_observed.second - l2_bias_.at(gnss_sat_id);

    double ionfree_phase = (pow(L1_frequency / L2_frequency, 2.0) * L1_lambda * L1_observed.first - L2_lambda * L2_observed.first) / (pow(L1_frequency / L2_frequency, 2.0) - 1);
    observed_values_.ionfree_carrier_phase.push_back(ionfree_phase);
  }
  return;
}

// この関数は外部で呼ばれて，内部が変更される．．．
void PBD_GnssObservation::UpdateInfoAfterObserved()
{
  // update observation state info
  for (int i = 0; i < num_of_gnss_satellites_; ++i)
  {
    // ここの操作はconstなのでできない．これはUpdate関数の初期部分に移行すればいい？
    info_.pre_observed_status.at(i) = info_.now_observed_status.at(i);
    info_.now_observed_status.at(i) = false;
  }
  info_.pre_observed_gnss_sat_id.clear();
  for (auto x : info_.now_observed_gnss_sat_id) info_.pre_observed_gnss_sat_id.push_back(x);

  return;
}

void PBD_GnssObservation::ClearPreValues(GnssObservedValues& values)
{
  values.observable_gnss_sat_id.clear();
  values.gnss_satellites_position.clear();
  values.gnss_clock.clear();
  values.L1_carrier_phase.clear();
  values.L1_pseudo_range.clear();
  values.L2_carrier_phase.clear();
  values.L2_pseudo_range.clear();
  values.ionfree_carrier_phase.clear();
  values.ionfree_pseudo_range.clear();
};

