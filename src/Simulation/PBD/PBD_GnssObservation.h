#ifndef PBD_GNSS_OBSERVATION_H_
#define PBD_GNSS_OBSERVATION_H_

#include <vector>
#include <utility>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include "Vector.hpp"
#include "GnssSatellites.h"
#include "./Orbit/Orbit.h"


// このクラスは受信機クラスか何かに拡張する．
struct GnssObservedValues
{
  std::vector<int> observable_gnss_sat_id;
  std::vector<libra::Vector<3>> gnss_satellites_position;
  std::vector<double> gnss_clock;

  std::vector<double> L1_pseudo_range; //[m]
  std::vector<double> L2_pseudo_range; //[m]
  std::vector<std::pair<double, double>> L1_carrier_phase; //[位相, bias(整数, dtype = double)]
  std::vector<std::pair<double, double>> L2_carrier_phase; //[位相, bias(整数, dtype = double)]

  std::vector<double> ionfree_pseudo_range;
  std::vector<double> ionfree_carrier_phase; //今入っているのは加工済み

  bool check_normal();
};

struct GnssObserveInfo
{
  // {gnss_sat_id: true/false, ...}
  std::vector<bool> pre_observed_status{};
  std::vector<bool> now_observed_status{};
  // {index: gnss_sat_id} index means the position in the state variables
  std::vector<int> pre_observed_gnss_sat_id{};
  std::vector<int> now_observed_gnss_sat_id{};
  // {index: info}
  // 以下はモデルで計算するものなのでここではない．
  // vector<double> geometric_range{};
  // vector<double> pseudo_range_model{};
  // vector<double> carrier_phase_range_model{};
};


// 普通にアクセスするよりGetter,使ってした方がいい．
class PBD_GnssObservation
{
public:
  PBD_GnssObservation(const Orbit& orbit, const GnssSatellites& gnss_satellites);
  ~PBD_GnssObservation();

  void Update();
  void UpdateGnssObservation();
  void CalcIonfreeObservation();
  void UpdateInfoAfterObserved();

  GnssObservedValues true_values_; // trueは要らんかも
  GnssObservedValues observed_values_;
  GnssObserveInfo info_;
  int num_of_gnss_satellites_; // というかこれはここに要らんのでは？
  std::vector<double> l1_bias_{};
  std::vector<double> l2_bias_{};

  // receiver clock biasの真値[m]
  double receiver_clock_bias_;

  //マスク角 [rad] <- これは衛星ごとに異なることが想定されるのでiniファイルとかで指定すべきでは？
  const double mask_angle = 10.0 / 180.0 * M_PI;

private:
  const Orbit& orbit_;
  const GnssSatellites& gnss_satellites_;
  //アンテナの中心の向きが、常に反地球方向を向いているとして、適当にマスク角を取って、その中にいるとする
  bool CheckCanSeeSatellite(const libra::Vector<3> satellite_position, const libra::Vector<3> gnss_position) const;
  void ClearPreValues(GnssObservedValues& values);

  // std::random_device seed_gen;
  std::mt19937 mt;
};


#endif
