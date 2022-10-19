#ifndef __PBD_dgps_H__
#define __PBD_dgps_H__

#include "Dense" // Eigen

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
#include "PBD_const.h"
#include <Environment/Global/CelestialRotation.h>

class PBD_dgps
{
public:
  PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Dynamics& main_dynamics, const Dynamics& target_dynamics, PBD_GnssObservation& main_observation, PBD_GnssObservation& target_observation, PBD_GeoPotential* geop); // OrbitとGnssObservation同時に取得したい．
  ~PBD_dgps();
  void Update(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, PBD_GnssObservation& main_observation, PBD_GnssObservation& target_observation, const CelestialRotation earth_rotation);//, const Orbit& main_orbit, const Orbit& target_orbit);
  void OrbitPropagation();
  //void calculate_difference_observation(GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, const int sat_id, Eigen::MatrixXd& pre_M);
  void KalmanFilter();
  void RemoveRows(Eigen::MatrixXd& matrix, unsigned int begin_row, unsigned int end_row);
  void RemoveColumns(Eigen::MatrixXd& matrix, unsigned int begin_col, unsigned int end_col);

private:
  PBD_GeoPotential* geo_potential_;

  // main satellite
  const Dynamics& main_dynamics_;
  // target satellite
  const Dynamics& target_dynamics_;

  // FIXME:データ構造は要検討．この構造にしていいのか？
  struct Ambiguity
  {
    std::vector<double> N; // [cycle] ambiguity
    std::vector<int> gnss_sat_id; // 対応するGNSS ID
    std::vector<bool> is_fixed; //不定性解除に成功したかどうか
  };

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
  };

  EstimatedVariables x_est_main;
  EstimatedVariables x_est_target;

  Eigen::VectorXd true_N_main;
  Eigen::VectorXd true_N_target; // [m]

  // std::map<const int, int> common_index_dict;

  // Eigen::VectorXd x_true; // 状態量真値 for log
  // シンプルにここが参照になっていないからか．
  std::vector<EstimatedVariables> x_est_; // 状態量ベクトル

  std::vector<libra::Vector<3>> antenna_pos_b_;
  std::vector<PBD_GnssObservation> gnss_observations_{};

  struct GnssObserveModel
  {
    std::vector<double> geometric_range;
    std::vector<double> pseudo_range_model;
    std::vector<double> carrier_phase_range_model;
  };
  std::vector<GnssObserveModel> gnss_observed_models_; // = { &main_model, &target_model };

  std::vector<int> visible_gnss_nums_{0, 0, 0}; // main, target, common
  std::vector<int> pre_visible_gnss_nums_{0, 0}; // main, target


  //初期化をもっとスマートにできるように考える
  //ここら辺も構造体にまとめるか．
  std::vector<bool> common_observed_status{};
  std::vector<int> common_observed_gnss_sat_id{};

  // air drag balistic coeficient
  const double Cd = 2.928e-14; // 高度に応じて変更したいが，高度変化ないから一旦，一定で行く．

  Eigen::MatrixXd P_;
  Eigen::MatrixXd Phi_;
  Eigen::MatrixXd Q_;
  Eigen::MatrixXd R_;

  double step_time;
  double observe_step_time = 10.0;
  double log_step_time = 1.0;

  libra::Matrix<3, 3> trans_eci_to_ecef_;

  // log用
  std::ofstream ofs;

  int num_of_gnss_satellites_;

  // receiver clock biasの真値[m] <- これはログ用にしか使ってないからイランかも
  const double& receiver_clock_bias_main_;
  const double& receiver_clock_bias_target_;

  std::random_device seed_gen;
  std::mt19937 mt;


  void InitAmbiguity(EstimatedVariables& x_est);
  // sat_idだけ指定する形でもいいかも
  void RK4(Eigen::Vector3d& position, Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration, Eigen::Vector3d& acc_dist, Eigen::MatrixXd& Phi);
  Eigen::Vector3d PositionDifferential(const Eigen::Vector3d& velocity) const;
  Eigen::Vector3d VelocityDifferential(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration, Eigen::Vector3d& acc_dist) const;
  void AddGeoPotentialDisturbance(const Eigen::Vector3d& position, Eigen::Vector3d& acc_dist) const;
  Eigen::MatrixXd UpdateP(void);
  Eigen::MatrixXd CalculateJacobian(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const;
  // Eigen::MatrixXd CalculateA(const Eigen::Vector3d& position_main, const Eigen::Vector3d& velocity_main, const Eigen::Vector3d& position_target, const Eigen::Vector3d& velocity_target);
  Eigen::Matrix3d TransRTN2ECI(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const;
  Eigen::MatrixXd CalculateQ_at(void); // 名前は要検討．
  void InitializeQ();
  Eigen::MatrixXd CalculatePhi_a(const double dt);
  Eigen::MatrixXd CalculateK(Eigen::MatrixXd H, Eigen::MatrixXd S);
  void ResizeS(Eigen::MatrixXd& S, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c);
  void ResizeMHt(Eigen::MatrixXd& MHt, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c);
  void UpdateTrueAmbiguity(std::vector<std::vector<double>> N, const int gnss_sat_id, const double lambda);
  void UpdateObservationsGRAPHIC(const int sat_id, EstimatedVariables& x_est, const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv, const int N_offset);
  void UpdateObservationsSDCP(const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv);
  void UpdateObservations(Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv);
  void FindCommonObservedGnss(const std::pair<int, int> sat_id_pair);
  void UpdateBiasForm(const int sat_id, EstimatedVariables& x_est, Eigen::MatrixXd& P, Eigen::MatrixXd& Q);
  void AllocateToCh(const int gnss_sat_id, std::map<const int, int>& observing_ch, std::vector<int>& free_ch);
  void RemoveFromCh(const int gnss_sat_id, std::map<const int, int>& observing_ch, std::vector<int>& free_ch);
  Eigen::VectorXd ConvReceivePosToCenterOfMass(Eigen::VectorXd x_state);
  Eigen::VectorXd ConvCenterOfMassToReceivePos(Eigen::VectorXd x_state);
  void InitLogTable(void);
  void PBD_dgps::InitializePhi(void);
  void ClearGnssObserveModels(GnssObserveModel& observed_model);
  void TransECI2RTN_P(Eigen::MatrixXd& P, Eigen::Matrix3d trans_eci_to_rtn);

  // Eigen::VectorXd CalculateSingleDifference(const Eigen::VectorXd& main_observation, const Eigen::VectorXd& target_observation) const;
  // 一旦singleだけにする
  //double calculate_double_difference(const Eigen::VectorXd& main_observation, const Eigen::VectorXd& target_observation) const;

  template <typename T> bool CheckVectorEqual(const std::vector<T>& a, const std::vector<T>& b);

  // void MakeDoubleDifference();
  int SelectBaseGnssSatellite(Eigen::VectorXd N, Eigen::MatrixXd P_N);
  void DynamicNoiseScaling(Eigen::MatrixXd Q_dash, Eigen::MatrixXd H);
};
#endif