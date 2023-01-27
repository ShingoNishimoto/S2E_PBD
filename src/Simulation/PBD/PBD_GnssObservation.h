#pragma once

#include <vector>
#include <utility>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include "Vector.hpp"
#include "GnssSatellites.h"
#include "./Orbit/Orbit.h"
#include "../Spacecraft/PBD_Components.h"

#define RANGE_OBSERVE_DEBUG

// ホンマは周波数帯とかコードの種類を抽象化したクラスを作ったほうがいいのかなあ．．
// このクラスは受信機クラスか何かに拡張する．
struct GnssObservedValues
{
  std::vector<int> observable_gnss_sat_id;
  std::vector<Vector<3>> gnss_satellites_position;
  std::vector<double> gnss_clock;

  std::vector<double> L1_pseudo_range; //[m]
  std::vector<double> L2_pseudo_range; //[m]
  std::vector<std::pair<double, double>> L1_carrier_phase; //[位相, N(整数, dtype = double)]
  std::vector<std::pair<double, double>> L2_carrier_phase; //[位相, N(整数, dtype = double)]

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
  std::vector<int> now_observed_gnss_sat_id{}; // こいつらのgetter作ってもいいな．．．
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
  // FIXME: constになおす．
  PBD_GnssObservation(PBD_GNSSReceiver* gnss_receiver, const GnssSatellites& gnss_satellites);
  ~PBD_GnssObservation();

  void Update(void);
  void UpdateGnssObservation();
  void UpdateInfoAfterObserved();

  double CalculatePseudoRange(const int gnss_sat_id, const Vector<3> sat_position, const double sat_clock) const;
  double CalculateCarrierPhase(const int gnss_sat_id, const Vector<3> sat_position, const double sat_clock, const double integer_bias, const double lambda, const double pcc) const;
  double CalculateGeometricRange(const int gnss_sat_id, const Vector<3> rec_position) const;
  double CalculateIonDelay(const int gnss_id, const Vector<3> rec_position, const double frequency) const; // GnssSatelliteとflag以外はIFをそろえている．

  // ここら辺を介す構成はやめたい．
  inline const Vector<3> GetAntennaAlignmentError(void) { return receiver_->GetAlignmentError(); }
  inline const Vector<3> GetAntennaPosition (void) { return receiver_->GetAntennaPositionBody(); }
  // inline GnssObserveInfo GetObserveInfo(void) const {return info_;}
  inline const int GetNowVisibleGnssNum(void) const {return info_.now_observed_gnss_sat_id.size();}
  inline const int GetPreVisibleGnssNum(void) const {return info_.pre_observed_gnss_sat_id.size();}
  inline const double GetGnssElevationDeg(const int ch) const {return receiver_->GetGnssInfo(ch).latitude * libra::rad_to_deg;}
  inline const double GetGnssAzimuthDeg(const int ch) const {return receiver_->GetGnssInfo(ch).longitude * libra::rad_to_deg;}
  const Vector<3> GetGnssDirection_c(const int ch) const;
  inline const Vector<3> GetGnssDirection_b(const int ch) const { return receiver_->TransCompoToBody(GetGnssDirection_c(ch)); }
  inline const Vector<3> GetGnssDirection_i(const int ch) { return receiver_->TransCompoToEci(GetGnssDirection_c(ch)); }
  inline const PBD_GNSSReceiver* GetReceiver(void) const { return receiver_; }

  GnssObservedValues true_values_; // trueは要らんかも
  GnssObservedValues observed_values_;
  GnssObserveInfo info_;
  int num_of_gnss_satellites_;
  std::vector<double> l1_bias_; // 初期捕捉時に確定するバイアス成分（捕捉している間は固定）
  std::vector<double> l2_bias_; // 初期捕捉時に確定するバイアス成分（捕捉している間は固定）

  //マスク角 [rad] <- これは衛星ごとに異なることが想定されるのでiniファイルとかで指定すべきでは？
  const double mask_angle = 10.0 / 180.0 * M_PI;

private:
  // const Orbit& orbit_;
  PBD_GNSSReceiver* receiver_;
  const GnssSatellites& gnss_satellites_;
  void ClearPreValues(GnssObservedValues& values);
  void ProcessGnssObservations(void);
  const double GetRangeComToArpTrue(const int ch) const;
  const double GetRangeComToArpDesign(const int ch) const;

  // std::random_device seed_gen;
  std::mt19937 mt;
};
