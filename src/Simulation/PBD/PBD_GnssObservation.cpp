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

PBD_GnssObservation::PBD_GnssObservation(const Orbit& orbit, const GnssSatellites& gnss_satellites) : orbit_(orbit), gnss_satellites_(gnss_satellites), num_of_gnss_satellites_(gnss_satellites_.GetNumOfSatellites())
{
  // 初期化
  gnss_observe_info_.pre_observed_status.assign(num_of_gnss_satellites_, false);
  gnss_observe_info_.now_observed_status.assign(num_of_gnss_satellites_, false);
}

PBD_GnssObservation::~PBD_GnssObservation() {}

void PBD_GnssObservation::Update(double sat_clock_true)
{

}

void PBD_GnssObservation::GetGnssPositionObservation(const double sat_clock_true)
{
  //推定値の計算
  num_of_gnss_satellites_ = gnss_satellites_.GetNumOfSatellites(); // 更新
  libra::Vector<3> sat_position_i = orbit_.GetSatPosition_i(); // ECI

  // この観測している衛星群の情報更新とかはここでした方がいいな．．．
  gnss_observe_info_.now_observed_gnss_sat_id.clear(); //クリア

  for (int gnss_sat_id = 0; gnss_sat_id < num_of_gnss_satellites_; ++gnss_sat_id) {
    //if(gnss_sat_id == 7 || gnss_sat_id == 23 || gnss_sat_id == 31) continue; 　←この衛星たちの軌道情報が悪いからこうしていたのか？
    if (!gnss_satellites_.GetWhetherValid(gnss_sat_id)) continue;
    libra::Vector<3> gnss_position = gnss_satellites_.Get_true_info().GetSatellitePositionEci(gnss_sat_id);
    bool see_flag = CheckCanSeeSatellite(sat_position_i, gnss_position);
    // init
    // main_index_dict.insert(std::make_pair(gnss_sat_id, -1));
    // common_index_dict.insert(std::make_pair(gnss_sat_id, -1));

    if (!see_flag)
    {
      // pre_main_observing_ch = now_main_observing_ch;
      // RemoveFromCh(gnss_sat_id, now_main_observing_ch, main_free_ch);
      continue;
    }
    observe_info.at(sat_id).now_observed_status.at(gnss_sat_id) = true;
    observe_info.at(sat_id).now_observed_gnss_sat_id.push_back(gnss_sat_id);
    int observed_gnss_index = observe_info.at(sat_id).now_observed_gnss_sat_id.size() - 1;
    if (sat_id == 0)
    {
      // main_index_dict.at(gnss_sat_id) = observed_gnss_index;
      // pre_main_observing_ch = now_main_observing_ch;
      //AllocateToCh(gnss_sat_id, now_main_observing_ch, main_free_ch);
    }
    double gnss_clock = gnss_satellites_.Get_true_info().GetSatelliteClock(gnss_sat_id); // これはclock bias
    double l1_pseudo_range = gnss_satellites_.GetPseudoRangeECI(gnss_sat_id, sat_position_i, sat_clock_true, L1_frequency);
    double l2_pseudo_range = gnss_satellites_.GetPseudoRangeECI(gnss_sat_id, sat_position_i, sat_clock_true, L2_frequency);
    auto l1_carrier_phase = gnss_satellites_.GetCarrierPhaseECI(gnss_sat_id, sat_position_i, sat_clock_true, L1_frequency);
    auto l2_carrier_phase = gnss_satellites_.GetCarrierPhaseECI(gnss_sat_id, sat_position_i, sat_clock_true, L2_frequency);

    double ionfree_range = (pow(L1_frequency / L2_frequency, 2.0) * l1_pseudo_range - l2_pseudo_range) / (pow(L1_frequency / L2_frequency, 2.0) - 1);
    double ionfree_phase = (pow(L1_frequency / L2_frequency, 2.0) * L1_lambda * (l1_carrier_phase.first + l1_carrier_phase.second) - L2_lambda * (l2_carrier_phase.first + l2_carrier_phase.second)) / (pow(L1_frequency / L2_frequency, 2.0) - 1);

    gnss_true_values_.observable_gnss_sat_id.push_back(gnss_sat_id);
    gnss_true_values_.gnss_satellites_position.push_back(gnss_position);
    gnss_true_values_.gnss_clock.push_back(gnss_clock);
    gnss_true_values_.L1_pseudo_range.push_back(l1_pseudo_range);
    gnss_true_values_.L2_pseudo_range.push_back(l2_pseudo_range);
    gnss_true_values_.L1_carrier_phase.push_back(l1_carrier_phase);
    gnss_true_values_.L2_carrier_phase.push_back(l2_carrier_phase);

    gnss_true_values_.ionfree_pseudo_range.push_back(ionfree_range);
    gnss_true_values_.ionfree_carrier_phase.push_back(ionfree_phase);

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

    gnss_observed_values_.observable_gnss_sat_id.push_back(gnss_sat_id);
    gnss_observed_values_.gnss_satellites_position.push_back(gnss_position);
    gnss_observed_values_.gnss_clock.push_back(gnss_clock);
    gnss_observed_values_.L1_pseudo_range.push_back(l1_pseudo_range);
    gnss_observed_values_.L2_pseudo_range.push_back(l2_pseudo_range);
    gnss_observed_values_.L1_carrier_phase.push_back(l1_carrier_phase);
    gnss_observed_values_.L2_carrier_phase.push_back(l2_carrier_phase);

    gnss_observed_values_.ionfree_pseudo_range.push_back(ionfree_range);
    gnss_observed_values_.ionfree_carrier_phase.push_back(ionfree_phase);
  }
}


bool PBD_GnssObservation::CheckCanSeeSatellite(const libra::Vector<3> satellite_position, const libra::Vector<3> gnss_position) const
{
  double angle_rad = angle(satellite_position, gnss_position - satellite_position);
  if (angle_rad < M_PI / 2.0 - mask_angle) return true;
  else return false;
}
