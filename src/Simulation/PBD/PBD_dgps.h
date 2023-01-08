#ifndef __PBD_dgps_H__
#define __PBD_dgps_H__

#include <Dense> // Eigen

#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <map>

#include "Vector.hpp"
#include "SimTime.h"
#include "GnssSatellites.h"
#include "./Orbit/Orbit.h"
#include "PBD_GnssObservation.h"
#include "PBD_GeoPotential.hpp"
#include "../Spacecraft/PBD_Sat.h"
#include "PBD_const.h"
#include <Environment/Global/CelestialRotation.h>
#include "PCCEstimation.hpp"
#include "PBD_Lambda.h"

class PBD_dgps
{
public:
  PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const std::vector<PBD_Sat*> spacecrafts, PBD_GeoPotential* geop); // OrbitとGnssObservation同時に取得したい．
  ~PBD_dgps();
  void Update(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, PBD_GnssObservation& main_observation, PBD_GnssObservation& target_observation, const CelestialRotation earth_rotation);//, const Orbit& main_orbit, const Orbit& target_orbit);

  void OrbitPropagation();
  void KalmanFilter();

  // publicにする意味ないな．
  void RemoveRows(Eigen::MatrixXd& matrix, unsigned int begin_row, unsigned int end_row);
  void RemoveColumns(Eigen::MatrixXd& matrix, unsigned int begin_col, unsigned int end_col);

private:
  PBD_GeoPotential* geo_potential_;

  struct EstimatedVariables
  {
    //[x[m], y[m], z[m]]
    Eigen::Vector3d position;
    // [cdt[m]]
    Eigen::VectorXd clock; // 複数GNSSなどになった時はvectorになりうる．
    //[vx[m/s], vy[m/s], vz[m/s]]
    Eigen::Vector3d velocity;
    //[ar[nm/s^2], at[nm/s^2], an[nm/s^2]]
    Eigen::Vector3d acceleration; // 経験加速度 RTN
    Eigen::Vector3d acc_dist; // 外乱による加速度（ログ用）ECI
    Ambiguity ambiguity;
    const double lambda = L1_lambda; // wave length [m]

    // Antenna parameters
    PhaseCenterCorrection* pcc;
  };

  EstimatedVariables x_est_main;
  EstimatedVariables x_est_target;
  PCCEstimation pcc_estimate_; // 相対的な量を推定するのでインスタンスとしては2衛星で一つ

  Eigen::VectorXd true_N_main;
  Eigen::VectorXd true_N_target; // [m]

  // std::map<const int, int> common_index_dict;

  // Eigen::VectorXd x_true; // 状態量真値 for log
  // シンプルにここが参照になっていないからか．
  std::vector<EstimatedVariables*> x_est_; // 状態量ベクトル

  std::vector<libra::Vector<3>> antenna_pos_b_;
  std::vector<PBD_GnssObservation> gnss_observations_{};

  struct GnssObserveModel
  {
    std::vector<double> geometric_range;
    std::vector<double> pseudo_range_model;
    std::vector<double> carrier_phase_model;
  };
  std::vector<GnssObserveModel> gnss_observed_models_; // = { &main_model, &target_model };

  // 参照したい情報をまとめたもの．この配列を持つことで参照して保持する．
  struct SatelliteInfo
  {
    const Dynamics& dynamics;
    EstimatedVariables& x_est;
    Eigen::VectorXd true_N;
    const PBD_Components* components;
  };
  std::vector<SatelliteInfo> sat_info_;

  std::vector<int> visible_gnss_nums_{0, 0, 0}; // main, target, common
  std::vector<int> pre_visible_gnss_nums_{0, 0, 0}; // main, target

  //初期化をもっとスマートにできるように考える
  //ここら辺も構造体にまとめるか．
  std::vector<bool> common_observed_status{};
  std::vector<int> common_observed_gnss_sat_id_{};
  std::vector<int> pre_common_observed_gnss_sat_id_{};

  // air drag ballistic coefficient
  const double Cd = 2.928e-14 / 5; // 高度に応じて変更したいが，高度変化ないから一旦，一定で行く．これも推定した方がいい気はする．

  // for Kalman Filter
  Eigen::MatrixXd P_;
  Eigen::MatrixXd Phi_;
  Eigen::MatrixXd Q_;
  Eigen::MatrixXd R_;
  Eigen::VectorXd x_;
  Eigen::VectorXd z_;
  Eigen::VectorXd hx_;
  Eigen::MatrixXd H_;
  Eigen::VectorXd pre_Rv_;

  double step_time;
  double observe_step_time = 10.0;
  double log_step_time = 1.0;

  libra::Matrix<3, 3> trans_eci_to_ecef_;

  // log用
  std::ofstream ofs;

  int num_of_gnss_satellites_;
  int num_observables_ = 0;
  int num_main_state_all_;
  int num_state_all_;

  std::random_device seed_gen;
  std::mt19937 mt;

  // 初期化関数
  void InitAmbiguity(EstimatedVariables& x_est);
  void InitializeQ();
  void InitLogTable(void);
  void PBD_dgps::InitializePhi(void);

  // Filter関連
  // sat_idだけ指定する形でもいいかも
  void RK4(Eigen::Vector3d& position, Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration, Eigen::Vector3d& acc_dist, Eigen::MatrixXd& Phi);
  Eigen::Vector3d PositionDifferential(const Eigen::Vector3d& velocity) const;
  Eigen::Vector3d VelocityDifferential(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration, Eigen::Vector3d& acc_dist) const;
  void AddGeoPotentialDisturbance(const Eigen::Vector3d& position, Eigen::Vector3d& acc_dist) const;
  Eigen::MatrixXd UpdateP(void);
  Eigen::MatrixXd CalculateJacobian(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const;
  Eigen::MatrixXd CalculatePhi_a(const double dt);
  Eigen::MatrixXd CalculateQ_at(void); // 名前は要検討．
  Eigen::MatrixXd CalculateK(Eigen::MatrixXd H);
  void DynamicNoiseScaling(Eigen::MatrixXd Q_dash, Eigen::MatrixXd H);

  // 観測関連
  void UpdateObservationsGRAPHIC(const int sat_id, EstimatedVariables& x_est, const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& hx, Eigen::MatrixXd& H, Eigen::VectorXd& Rv, const int N_offset);
  void UpdateObservationsSDCP(const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& hx, Eigen::MatrixXd& H, Eigen::VectorXd& Rv);
  void UpdateObservations(Eigen::VectorXd& z, Eigen::VectorXd& hx, Eigen::MatrixXd& H, Eigen::VectorXd& Rv);
  void FindCommonObservedGnss(const std::pair<int, int> sat_id_pair);
  void UpdateBiasForm(const int sat_id, EstimatedVariables& x_est, Eigen::MatrixXd& P, Eigen::MatrixXd& Q);
  void ClearGnssObserveModels(GnssObserveModel& observed_model);
  void InitGnssObserveModels(GnssObserveModel& observed_model, const int gnss_num);
  Eigen::Vector3d ConvReceivePosToCenterOfMass(Eigen::Vector3d rec_pos, libra::Vector<3> antenna_pos_b, const Dynamics& dynamics);
  Eigen::Vector3d ConvCenterOfMassToReceivePos(Eigen::Vector3d pos, libra::Vector<3> antenna_pos_b, const Dynamics& dynamics);
  void AdjustReceiveCovariance(const std::vector<int>& now_gnss_sat_ids, const std::vector<int>& pre_gnss_sat_ids, const int gnss_sat_id, const int base_offset, const int pre_base_offset, const Eigen::VectorXd& pre_Rv);
  void UpdateNumOfState(PBD_GnssObservation main_observation, PBD_GnssObservation target_observation);
  const double DataEditing(const int sat_id, Eigen::VectorXd z);

// 整数不定性解除
  const bool IntegerAmbiguityResolution(const Eigen::VectorXd& x_update);

  // PCO, PCV関連
  const bool EstimateRelativePCC(const std::vector<double> sdcp_vec, const double elapsed_time);

  // 便利関数関連
  void TransECI2RTN_P(Eigen::MatrixXd& P, Eigen::Matrix3d trans_eci_to_rtn);
  Eigen::Matrix3d TransRTN2ECI(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const;
  void ResizeS(Eigen::MatrixXd& S, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c);
  void ResizeMHt(Eigen::MatrixXd& MHt, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c);
  template <typename T> bool CheckVectorEqual(const std::vector<T>& a, const std::vector<T>& b);
  void SetStateToVector(const int num_state_all, Eigen::VectorXd& x_state); // 代入用関数．
  void GetStateFromVector(const int num_main_state_all, const Eigen::VectorXd& x_state);

};
#endif