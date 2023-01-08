#include "PBD_GnssObservation.h"
#include "PBD_const.h"

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

PBD_GnssObservation::PBD_GnssObservation(PBD_GNSSReceiver* gnss_receiver, const GnssSatellites& gnss_satellites) : receiver_(gnss_receiver), gnss_satellites_(gnss_satellites)
{
  num_of_gnss_satellites_ = gnss_satellites_.GetNumOfSatellites();
  l1_bias_.assign(num_of_gnss_satellites_, 0.0);
  l2_bias_.assign(num_of_gnss_satellites_, 0.0);
  info_.pre_observed_status.assign(num_of_gnss_satellites_, false);
  info_.now_observed_status.assign(num_of_gnss_satellites_, false);
}

PBD_GnssObservation::~PBD_GnssObservation() {}

void PBD_GnssObservation::Update(void)
{
  UpdateGnssObservation();
  ProcessGnssObservations();
}

void PBD_GnssObservation::UpdateGnssObservation()
{
  //推定値の計算
  num_of_gnss_satellites_ = gnss_satellites_.GetNumOfSatellites(); // 時刻によって利用可能衛星数が変化するので更新
  // ここでMainRoutineを呼べばいいだけではあった．constのままでしたいな．

  const Vector<3> antenna_position_i = receiver_->GetAntennaPositionTrueECI(); // ECI

  info_.now_observed_status.assign(num_of_gnss_satellites_, false);
  info_.now_observed_gnss_sat_id.clear(); //クリア
  ClearPreValues(true_values_);
  ClearPreValues(observed_values_);

  const std::vector<GnssInfo> vec_gnss_info = receiver_->GetGnssInfoVec();

  // 受信機で受信できた衛星だけにアクセス．
  for (int ch = 0; ch < vec_gnss_info.size(); ch++)
  {
    const int gnss_sat_id = gnss_satellites_.GetIndexFromID(vec_gnss_info.at(ch).ID); // idとindexの定義を混同しないように整理する．

    info_.now_observed_status.at(gnss_sat_id) = true;
    info_.now_observed_gnss_sat_id.push_back(gnss_sat_id);

    PBD_GNSSReceiver::GnssReceiverObservations raw_observation = receiver_->GetRawObservations(ch);
    double l1_pseudo_range = raw_observation.l1_pseudo_range;
    double l2_pseudo_range = raw_observation.l2_pseudo_range;
    auto l1_carrier_phase  = raw_observation.l1_carrier_phase;
    auto l2_carrier_phase  = raw_observation.l2_carrier_phase;

// 受信機クラスでやってしまってもいい気はするが，受信機からは生のデータを得るというのが基本的な考えな気がするのでここで実施．
#ifdef RANGE_OBSERVE_DEBUG
    // ここで重心位置までの位置への変換は衛星側の処理として実施するので実際の値はわからず，設計値で変換する．
    const double range_com_to_arp = GetRangeComToArpDesign(ch);
    l1_pseudo_range += range_com_to_arp;
    l2_pseudo_range += range_com_to_arp;
    l1_carrier_phase.first += range_com_to_arp / L1_lambda;
    l2_carrier_phase.first += range_com_to_arp / L2_lambda;
#endif // RANGE_OBSERVE_DEBUG

    observed_values_.L1_pseudo_range.push_back(l1_pseudo_range);
    observed_values_.L2_pseudo_range.push_back(l2_pseudo_range);
    observed_values_.L1_carrier_phase.push_back(l1_carrier_phase);
    observed_values_.L2_carrier_phase.push_back(l2_carrier_phase);

    // 擬似距離の平均をとってクロック誤差を小さくするための処理をする．それをしないとそもそもなのでは？

    // true info
    Vector<3> gnss_position = gnss_satellites_.Get_true_info().GetSatellitePositionEci(gnss_sat_id);
    double gnss_clock = gnss_satellites_.Get_true_info().GetSatelliteClock(gnss_sat_id); // これはclock bias
    // estimateに使う方の情報
    gnss_position = gnss_satellites_.GetSatellitePositionEci(gnss_sat_id);
    gnss_clock = gnss_satellites_.GetSatelliteClock(gnss_sat_id);

    observed_values_.observable_gnss_sat_id.push_back(gnss_sat_id);
    observed_values_.gnss_satellites_position.push_back(gnss_position);
    observed_values_.gnss_clock.push_back(gnss_clock);
  }
}

void PBD_GnssObservation::ProcessGnssObservations(void)
{
  // ここの前半の処理はUpdateの方にまとめてもいいかも．
  int observed_gnss_index = 0;

  const std::vector<GnssInfo> vec_gnss_info = receiver_->GetGnssInfoVec();
  // 受信機で受信できた衛星だけにアクセス．
  for (int ch = 0; ch < vec_gnss_info.size(); ch++)
  {
    const int gnss_sat_id = gnss_satellites_.GetIndexFromID(vec_gnss_info.at(ch).ID); // idとindexの定義を混同しないように整理する．

    if (info_.pre_observed_status.at(gnss_sat_id) == true && info_.now_observed_status.at(gnss_sat_id) == false)
    {
      l1_bias_.at(gnss_sat_id) = 0.0;
      l2_bias_.at(gnss_sat_id) = 0.0;
    }
    else if (info_.pre_observed_status.at(gnss_sat_id) == false && info_.now_observed_status.at(gnss_sat_id) == true)
    {
      l1_bias_.at(gnss_sat_id) = observed_values_.L1_carrier_phase.at(observed_gnss_index).second;
      l2_bias_.at(gnss_sat_id) = observed_values_.L2_carrier_phase.at(observed_gnss_index).second;
    }
    if (info_.now_observed_status.at(gnss_sat_id)) ++observed_gnss_index;
  }

  for (int index = 0; index < info_.now_observed_gnss_sat_id.size(); index++)
  {
    int gnss_sat_id = info_.now_observed_gnss_sat_id.at(index);

    auto L1_observed = observed_values_.L1_carrier_phase.at(index);
    L1_observed.first += L1_observed.second - l1_bias_.at(gnss_sat_id); // 観測量には追尾分の波長変化も含める．
    observed_values_.L1_carrier_phase.at(index).first = L1_observed.first; // 位相観測量として更新

    auto L2_observed = observed_values_.L2_carrier_phase.at(index);
    L2_observed.first += L2_observed.second - l2_bias_.at(gnss_sat_id); // 観測量には追尾分の波長変化も含める．
    observed_values_.L2_carrier_phase.at(index).first = L2_observed.first; // 位相観測量として更新

    // double ionfree_phase = L2_lambda * (L1_frequency / L2_frequency * (L1_observed.first + l1_carrier_phase.second) - (l2_carrier_phase.first + l2_carrier_phase.second)) / (pow(L1_frequency / L2_frequency, 2.0) - 1);
  }
}

const Vector<3> PBD_GnssObservation::GetGnssDirection_c(const int ch) const
{
  const double azi_rad = receiver_->GetGnssInfo(ch).longitude;
  const double ele_rad = receiver_->GetGnssInfo(ch).latitude;
  Vector<3> e(0);
  e[0] = cos(ele_rad) * cos(azi_rad); e[1] = cos(ele_rad) * sin(azi_rad); e[2] = sin(ele_rad);
  return e;
}

// この関数は外部で呼ばれて，内部が変更される
void PBD_GnssObservation::UpdateInfoAfterObserved()
{
  // update observation state info
  info_.pre_observed_gnss_sat_id.clear();
  for (int index = 0; index < info_.now_observed_gnss_sat_id.size(); index++)
  {
    int gnss_sat_id = info_.now_observed_gnss_sat_id.at(index);
    info_.pre_observed_status.at(gnss_sat_id) = info_.now_observed_status.at(gnss_sat_id);
    info_.pre_observed_gnss_sat_id.push_back(gnss_sat_id);
  }
  // nowのクリアは次回更新の最初に行う．
  // info_.now_observed_status.assign(num_of_gnss_satellites_, false);
  // info_.now_observed_gnss_sat_id.clear();
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

double PBD_GnssObservation::CalculatePseudoRange(const int gnss_sat_id, const Vector<3> sat_position, const double sat_clock) const
{
  // この情報もテーブルとして持っておけば無駄がない．
  const int index = GetIndexOfStdVector<int>(info_.now_observed_gnss_sat_id, gnss_sat_id);
  double gnss_clock = observed_values_.gnss_clock.at(index);

  double range = 0.0;
  // 推定するときはここにアライメント誤差の推定量も混ぜる．
  Vector<3> receive_position = receiver_->GetCodeReceivePositionDesignECI(sat_position);
  range = CalculateGeometricRange(gnss_sat_id, receive_position);

  // clock offsetの分を追加
  range += sat_clock - gnss_clock; // 電離層はフリーにしている．

#ifdef RANGE_OBSERVE_DEBUG
  range += GetRangeComToArpDesign(index);
#endif // RANGE_OBSERVE_DEBUG

  return range;
}

double PBD_GnssObservation::CalculateCarrierPhase(const int gnss_sat_id, const Vector<3> sat_position,
        const double sat_clock, const double integer_bias, const double lambda, const double pcc) const
{
  const int index = GetIndexOfStdVector<int>(info_.now_observed_gnss_sat_id, gnss_sat_id);
  double gnss_clock = observed_values_.gnss_clock.at(index);

  double range = 0.0;
  Vector<3> receive_position = receiver_->GetPhaseReceivePositionDesignECI(sat_position);
  range = CalculateGeometricRange(gnss_sat_id, receive_position);

  range += sat_clock - gnss_clock;
  range += lambda * integer_bias; // ここも電離圏は入れてない．
  range += pcc;

#ifdef RANGE_OBSERVE_DEBUG
  range += GetRangeComToArpDesign(index);
#endif // RANGE_OBSERVE_DEBUG

  return range; // 位相観測量に変換（単位は[m]）
}

double PBD_GnssObservation::CalculateGeometricRange(const int gnss_sat_id, const Vector<3> rec_position) const
{
  const int index = GetIndexOfStdVector<int>(info_.now_observed_gnss_sat_id, gnss_sat_id);
  auto gnss_position = observed_values_.gnss_satellites_position.at(index);

  double range = 0.0;
  for (int i = 0; i < 3; ++i) {
    range += pow(rec_position[i] - gnss_position[i], 2.0);
  }
  range = sqrt(range);

  return range;
}

// 内容はほぼGnssSatelliteからのコピー
double PBD_GnssObservation::CalculateIonDelay(const int gnss_id, const Vector<3> rec_position, const double frequency) const
{
  // gnss_id is wrong or not validate
  if (gnss_id >= num_of_gnss_satellites_) return 0.0;

  const double Earth_hemisphere = 6378.1;  //[km]

  double altitude = 0.0;
  for (int i = 0; i < 3; ++i) altitude += pow(rec_position[i], 2.0);
  altitude = sqrt(altitude);
  altitude = altitude / 1000.0 - Earth_hemisphere;  //[m -> km]
  if (altitude >= 1000.0) return 0.0;               // there is no Ionosphere above 1000km

  Vector<3> gnss_position;
  gnss_position = gnss_satellites_.GetSatellitePositionEci(gnss_id);

  double angle_rad = angle(rec_position, gnss_position - rec_position);
  const double default_delay = 20.0;                                             //[m] default delay
  double delay = default_delay * (1000.0 - altitude) / 1000.0 / cos(angle_rad);  // set the maximum height as 1000.0. Divide by
                                                                                 // cos because the slope makes it longer.
  const double default_frequency = 1500.0;                                       //[MHz]
  // Ionospheric delay is inversely proportional to the square of the frequency
  delay *= pow(default_frequency / frequency, 2.0);

  return delay;
}

// 実際のアライメント誤差ありの情報
const double PBD_GnssObservation::GetRangeComToArpTrue(const int ch) const
{
  Vector<3> gnss_direction_b = GetGnssDirection_b(ch);
  const Vector<3> arp_true = receiver_->GetAntennaPositionBody() + receiver_->GetAlignmentError();
  const double range_com_to_arp = libra::inner_product(arp_true, gnss_direction_b);

  return range_com_to_arp;
}

const double PBD_GnssObservation::GetRangeComToArpDesign(const int ch) const
{
  Vector<3> gnss_direction_b = GetGnssDirection_b(ch);
  const double range_com_to_arp = libra::inner_product(receiver_->GetAntennaPositionBody(), gnss_direction_b);

  return range_com_to_arp;
}
