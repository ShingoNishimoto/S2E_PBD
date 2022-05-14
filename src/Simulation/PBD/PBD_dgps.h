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


class PBD_dgps
{
public:
  struct EstimatedVariables
  {
    //[x[m], y[m], z[m]]
    Eigen::Vector3d position;
    // [cdt[m]]
    Eigen::VectorXd clock; // 複数GNSSなどになった時はvectorになりうる．
    //[vx[m/s], vy[m/s], vz[m/s]]
    Eigen::Vector3d velocity;
    //[ax[nm/s^2], ay[nm/s^2], az[nm/s^2]]
    Eigen::Vector3d acceleration; // 経験加速度
    Eigen::VectorXd bias; // [m] ambiguityと書く方がいい？
  };

  PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Orbit& main_orbit, const Orbit& target_orbit, PBD_GnssObservation& main_observation, PBD_GnssObservation& target_observation); // OrbitとGnssObservation同時に取得したい．
  ~PBD_dgps();
  void Update(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_);//, const Orbit& main_orbit, const Orbit& target_orbit);
  void OrbitPropagation();
  void SetBiasToObservation(const int sat_id, EstimatedVariables& x_est, PBD_GnssObservation& gnss_observation);
  // そもそもこれをpublicにしている理由がないか．
  void UpdateBiasForm(const int sat_id, EstimatedVariables& x_est);
  //void calculate_difference_observation(GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, const int sat_id, Eigen::MatrixXd& pre_M);
  void KalmanFilter();
  void RemoveRows(Eigen::MatrixXd& matrix, unsigned int begin_row, unsigned int end_row);
  void RemoveColumns(Eigen::MatrixXd& matrix, unsigned int begin_col, unsigned int end_col);
  std::ofstream ofs;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  // 複数衛星に拡張するのならorbitsにした方がいい．
  // main satellite
  const Orbit& main_orbit_;
  // target satellite
  const Orbit& target_orbit_;

  EstimatedVariables x_est_main;
  EstimatedVariables x_est_target;

  Eigen::VectorXd true_bias_main;

  // std::map<const int, int> main_index_dict;

  // 辞書が欲しい commonとmainをつなげるために
  Eigen::VectorXd true_bias_target; // [m]

  // std::map<const int, int> common_index_dict;

  // Eigen::VectorXd x_true; // 状態量真値 for log
  // シンプルにここが参照になっていないからか．
  vector<EstimatedVariables> x_est{}; // 状態量ベクトル

  // gnss_sat_id <-> indexの変換が簡単にできるようにしたい．

  vector<PBD_GnssObservation*> gnss_observations_;

    struct GnssObserveModel
    {
      vector<double> geometric_range;
      vector<double> pseudo_range_model;
      vector<double> carrier_phase_range_model;
    };

    vector<GnssObserveModel*> gnss_observed_models_; // = { &main_model, &target_model }; //なぜここでの初期化が必要なのか?あと普通に参照型であれば無理なのはなぜ？

    //初期化をもっとスマートにできるように考える
    //ここら辺も構造体にまとめるか．
    vector<bool> common_observed_status{};
    vector<int> common_observed_gnss_sat_id{};
    
    // now, preが必要か？
    /*
    std::map<const int, int> pre_main_observing_ch; // <ch, gnss_sat_id>
    std::map<const int, int> now_main_observing_ch;
    std::map<const int, int> pre_common_observing_ch;
    std::map<const int, int> now_common_observing_ch;
    vector<int> main_free_ch{};
    vector<int> common_free_ch{};
    */

    // air drag balistic coeficient
    const double Cd = 2.928e-14; // 高度に応じて変更したいが，高度変化ないから一定でいいか．

    Eigen::MatrixXd M;

    const int num_of_single_status = 10;
    const int num_of_status = 20;
    const int num_of_gnss_channel = 12; // max receivable gnss number
    const int observation_dimension = 3*num_of_gnss_channel; // GRAPHIC*2 + SDCP
    const int state_dimension = num_of_status + 2 * num_of_gnss_channel;
    const int single_dimension = num_of_single_status + num_of_gnss_channel;

    double step_time;
    double observe_step_time = 10.0;
    double log_step_time = 1.0;
    vector<Eigen::Vector3d> RK4(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration);
    Eigen::Vector3d PositionDifferential(const Eigen::Vector3d& velocity) const;
    Eigen::Vector3d VelocityDifferential(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration) const;
    // for differential
    Eigen::MatrixXd UpdateM();
    // Eigen::MatrixXd CalculateA(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration) const;
    // for differential
    Eigen::MatrixXd CalculateA(const EstimatedVariables& x_est_main, const EstimatedVariables& x_est_target) const;
    Eigen::MatrixXd CalculateJacobian(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration) const;
    Eigen::Matrix3d TransRTN2ECI(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const;
    Eigen::MatrixXd CalculateQ();
    Eigen::MatrixXd CalculatePhi_a(const double dt);
    Eigen::MatrixXd CalculateK(Eigen::MatrixXd H, Eigen::MatrixXd S);
    void ResizeS(Eigen::MatrixXd& S, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c);
    void ResizeMHt(Eigen::MatrixXd& MHt, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c);
    void UpdateTrueBias(vector<vector<double>> bias, const int gnss_sat_id, const double lambda);
    void UpdateObservationsGRAPHIC(const int sat_id, EstimatedVariables& x_est, const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv);
    void UpdateObservationsSDCP(const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv);
    void FindCommonObservedGnss(const std::pair<int, int> sat_id_pair);
    void AllocateToCh(const int gnss_sat_id, std::map<const int, int>& observing_ch, vector<int>& free_ch);
    void RemoveFromCh(const int gnss_sat_id, std::map<const int, int>& observing_ch, vector<int>& free_ch);

    int num_of_gnss_satellites_;

    // receiver clock biasの真値[m] <- これはログ用にしか使ってないからイランかも
    const double& receiver_clock_bias_main_;
    const double& receiver_clock_bias_target_;

    std::random_device seed_gen;
    std::mt19937 mt;

    double CalculatePseudoRange(const EstimatedVariables& x_est, libra::Vector<3> gnss_position, double gnss_clock);
    double CalculateCarrierPhase(const EstimatedVariables& x_est, libra::Vector<3> gnss_position, double gnss_clock, double integer_bias, double lambda);
    double CalculateGeometricRange(const Eigen::Vector3d& sat_position, libra::Vector<3> gnss_position) const;
    Eigen::VectorXd CalculateSingleDifference(const Eigen::VectorXd& main_observation, const Eigen::VectorXd& target_observation) const;
    // 一旦singleだけにする
    //double calculate_double_difference(const Eigen::VectorXd& main_observation, const Eigen::VectorXd& target_observation) const;

    template <typename T> bool CheckVectorEqual(const vector<T>& a, const vector<T>& b);
};
#endif