#include <iomanip>
#include <set>
#include "PBD_dgps.h"
#include "PBD_Lambda.h"
#include <GeoPotential.h>
#include "../../Library/VectorTool.hpp"

// clock model
#define CLOCK_IS_WHITE_NOISE (1)
// #define CLOCK_IS_RANDOM_WALK (1)

// inverse calculation method
#define CHOLESKY (0)
#define QR (1)

// Kalman Filter method
#define AKF
#define AIR_DRAG_ON

#define REDUCED_DYNAMIC

#ifdef REDUCED_DYNAMIC
#define NUM_SINGLE_STATE (10)
#else
#define NUM_SINGLE_STATE (7)
#endif

#define NUM_SAT (2)
#define NUM_STATE (NUM_SINGLE_STATE * NUM_SAT)
#define NUM_GNSS_CH (12)
#define NUM_OBSERVABLES (NUM_GNSS_CH * 3)
#define NUM_SINGLE_STATE_ALL (NUM_SINGLE_STATE + NUM_GNSS_CH)
#define NUM_STATE_ALL (NUM_SINGLE_STATE_ALL * NUM_SAT)

// #define GEOP_DEBUG
// #define N_DEBUG
// Nの真値をしっかり実装しないとダメかも．
// #define TIME_UPDATE_DEBUG

#undef cross

static const int conv_index_from_gnss_sat_id(std::vector<int> observed_gnss_sat_id, const int gnss_sat_id);
static const double PBD_DGPS_kConvNm2m = 1e-9;
static const int precision = 10; // position

// template使っていい感じにしたい．
// template <typename T> static void LogOutput(const )
static void LogOutput_(std::ofstream& ofs, const Eigen::MatrixXd& M, const int size, const int max_size);

// outputを変えるときは"result.csv"を変更する．せめてパスは変えたい．
PBD_dgps::PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Dynamics& main_dynamics, const Dynamics& target_dynamics, PBD_GnssObservation& main_observation, PBD_GnssObservation& target_observation, PBD_GeoPotential* geop) :mt(42), step_time(sim_time_.GetStepSec()), ofs("result_new.csv"), num_of_gnss_satellites_(gnss_satellites_.GetNumOfSatellites()), main_dynamics_(main_dynamics), target_dynamics_(target_dynamics), gnss_observations_({ main_observation, target_observation }), geo_potential_(geop)
{
  //初期化
  x_est_main.position = Eigen::VectorXd::Zero(3);
  libra::Vector<3> position_main = main_dynamics_.GetPosition_i();
  for (int i = 0; i < 3; ++i) x_est_main.position(i) = position_main[i];
  x_est_main.clock = Eigen::VectorXd::Zero(1);
  libra::Vector<3> velocity_main = main_dynamics_.GetOrbit().GetSatVelocity_i();
  for (int i = 0; i < 3; ++i) x_est_main.velocity(i) = velocity_main[i];
  x_est_main.acceleration = Eigen::VectorXd::Zero(3);
  x_est_main.acc_dist = Eigen::VectorXd::Zero(3);
  InitAmbiguity(x_est_main);

  x_est_target.position = Eigen::VectorXd::Zero(3);
  libra::Vector<3> position_target = target_dynamics_.GetPosition_i();
  for (int i = 0; i < 3; ++i) x_est_target.position(i) = position_target[i];
  x_est_target.clock = Eigen::VectorXd::Zero(1);
  libra::Vector<3> velocity_target = target_dynamics_.GetOrbit().GetSatVelocity_i();
  for (int i = 0; i < 3; ++i) x_est_target.velocity(i) = velocity_target[i];
  x_est_target.acceleration = Eigen::VectorXd::Zero(3);
  x_est_target.acc_dist = Eigen::VectorXd::Zero(3);
  InitAmbiguity(x_est_target);

  // x_est_.push_back(x_est_main);
  // x_est_.push_back(x_est_target);
  // true_N_main = Eigen::VectorXd::Zero(NUM_GNSS_CH);
  // true_N_target = Eigen::VectorXd::Zero(NUM_GNSS_CH);
  sat_info_.push_back({main_dynamics, x_est_main, Eigen::VectorXd::Zero(NUM_GNSS_CH)});
  sat_info_.push_back({target_dynamics, x_est_target, Eigen::VectorXd::Zero(NUM_GNSS_CH)});

  GnssObserveModel main_model{};
  GnssObserveModel target_model{};
  gnss_observed_models_.push_back(main_model);
  gnss_observed_models_.push_back(target_model);

  antenna_pos_b_.push_back(main_observation.GetAntennaPosition());
  antenna_pos_b_.push_back(target_observation.GetAntennaPosition());

  // STMの初期化
  InitializePhi();

  // 初期分散
  std::normal_distribution<> position_dist(0.0,sigma_r_ini);
  std::normal_distribution<> velocity_dist(0.0, sigma_v_ini);
  std::normal_distribution<> acc_r_dist(0.0, sigma_acc_r_ini);
  std::normal_distribution<> acc_t_dist(0.0, sigma_acc_t_ini);
  std::normal_distribution<> acc_n_dist(0.0, sigma_acc_n_ini);
  std::normal_distribution<> N_dist(0.0, sigma_N_ini);

  Eigen::VectorXd V = Eigen::VectorXd::Zero(NUM_STATE);

  // 状態量を減らした部分を実装.
  // A-priori
  for(int i = 0; i < 3; ++i) V(i) = pow(sigma_r_ini, 2.0); // position
  V(3) = pow(sigma_cdt_ini, 2.0); // clock
  for(int i = 0; i < 3; ++i) V(4 + i) = pow(sigma_v_ini, 2.0); // velocity
  // acceleration
#ifdef REDUCED_DYNAMIC
  V(7) = pow(sigma_acc_r_ini, 2.0);
  V(8) = pow(sigma_acc_t_ini, 2.0);
  V(9) = pow(sigma_acc_n_ini, 2.0);
#endif // REDUCED_DYNAMIC

  // 以下がtarget
  for (int i = 0; i < 3; ++i) V(NUM_SINGLE_STATE + i) = pow(sigma_r_ini, 2.0);
  V(NUM_SINGLE_STATE + 3) = pow(sigma_cdt_ini, 2.0);
  for (int i = 0; i < 3; ++i) V(NUM_SINGLE_STATE + 4 + i) = pow(sigma_v_ini, 2.0);
#ifdef REDUCED_DYNAMIC
  V(NUM_SINGLE_STATE + 7) = pow(sigma_acc_r_ini, 2.0);
  V(NUM_SINGLE_STATE + 8) = pow(sigma_acc_t_ini, 2.0);
  V(NUM_SINGLE_STATE + 9) = pow(sigma_acc_n_ini, 2.0);
#endif // REDUCED_DYNAMIC

  P_ = V.asDiagonal(); // 誤差分散共分散行列P 初めはN無し

  // Process noise Q_
  InitializeQ();
  // visible_gnss_nums_ = {0, 0, 0};

  // Measurement noise R_
  Eigen::VectorXd Rv = Eigen::VectorXd::Zero(NUM_OBSERVABLES);
  Rv.topRows(2*NUM_GNSS_CH) = (pow(pseudo_sigma / 2.0, 2.0) + pow(carrier_sigma / 2.0, 2.0))* Eigen::VectorXd::Ones(2*NUM_GNSS_CH); // GRAPHIC
  Rv.bottomRows(NUM_GNSS_CH) = 2.0*pow(carrier_sigma, 2.0) * Eigen::VectorXd::Ones(NUM_GNSS_CH); // SDCP
  R_ = Rv.asDiagonal();

  // 初期位置はガウシアンからサンプリング．mtは乱数のシード
  for(int i = 0; i < 3; ++i) x_est_main.position(i) += position_dist(mt);
  for (int i = 0; i < 3; ++i) x_est_main.velocity(i) += velocity_dist(mt);

  for(int i = 0; i < 3; ++i) x_est_target.position(i) += position_dist(mt);
  for (int i = 0; i < 3; ++i) x_est_target.velocity(i) += velocity_dist(mt);

  common_observed_status.assign(num_of_gnss_satellites_, false);


  std::ofstream ofs_ini_txt("readme_new.txt");
  ofs_ini_txt << "initial position dist: " << sigma_r_ini << std::endl;
  ofs_ini_txt << "initial velocity dist: " << sigma_v_ini << std::endl;
  ofs_ini_txt << "initial acceleration dist: " << sigma_acc_r_ini << std::endl;
  ofs_ini_txt << "initial clock dist: " << sigma_cdt_ini << std::endl;
  ofs_ini_txt << "initial ambiguity dist[cycle]: " << sigma_N_ini << std::endl;
  ofs_ini_txt << "pseudo dist: " << pseudo_sigma << std::endl;
  ofs_ini_txt << "carrier dist: " << carrier_sigma << std::endl;
  ofs_ini_txt << "clock dist: " << clock_sigma << std::endl;
  ofs_ini_txt << "process noise of position: " << sigma_r_process << std::endl;
  ofs_ini_txt << "process noise of velocity: " << sigma_v_process << std::endl;
  ofs_ini_txt << "process noise of radial acceleration: " << sigma_acc_r_process << std::endl;
  ofs_ini_txt << "process noise of tangential acceleration: " << sigma_acc_t_process << std::endl;
  ofs_ini_txt << "process noise of north acceleration: " << sigma_acc_n_process << std::endl;
  ofs_ini_txt << "process noise of clock: " << sigma_cdt_process << std::endl;
  ofs_ini_txt << "process noise of ambiguity[cycle]: " << sigma_N_process << std::endl;
  ofs_ini_txt << "time const. acceleration: " << tau_a << std::endl;
  ofs_ini_txt << "time const. clock: " << tau_cdt << std::endl;
  ofs_ini_txt << "mask angle: " << gnss_observations_.at(0).mask_angle << std::endl; // FIXME
  ofs_ini_txt << "num of status: " << NUM_STATE << std::endl;
  ofs_ini_txt << "observe step time: " << observe_step_time << std::endl;
  ofs_ini_txt << "log step time: " << log_step_time << std::endl;
  libra::Vector<3> alignment_err_main = main_observation.GetAntennaAlignmentError();
  libra::Vector<3> alignment_err_target = target_observation.GetAntennaAlignmentError();
  for (uint8_t i = 0; i < 3; i++)
    ofs_ini_txt << "main antenna alignment error[m]: " << alignment_err_main[i] << std::endl;
  for (uint8_t i = 0; i < 3; i++)
    ofs_ini_txt << "target antenna alignment error[m]: " << alignment_err_target[i] << std::endl;
}

PBD_dgps::~PBD_dgps(){}

void PBD_dgps::InitAmbiguity(EstimatedVariables& x_est)
{
  x_est.ambiguity.N.assign(NUM_GNSS_CH, 0);
  x_est.ambiguity.gnss_sat_id.assign(NUM_GNSS_CH, num_of_gnss_satellites_);
  x_est.ambiguity.is_fixed.assign(NUM_GNSS_CH, false);
}

// ここが他と同じ時刻系を使ってないのが原因な気がしてきた．FIXME！
void PBD_dgps::Update(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, PBD_GnssObservation& main_observation_, PBD_GnssObservation& target_observation_, const CelestialRotation earth_rotation)
{
  trans_eci_to_ecef_ = earth_rotation.GetDCMJ2000toXCXF();

  double elapsed_time = sim_time_.GetElapsedSec();
  double tmp = floor(elapsed_time/observe_step_time + 1e-4); //1e-4は数値誤差
  double tmp_log = floor(elapsed_time/log_step_time + 1e-4);

  //まず更新
  OrbitPropagation();

  //観測時間にピッタリ
  if(abs(elapsed_time - tmp*observe_step_time) < step_time/2.0){
#ifndef TIME_UPDATE_DEBUG
    P_ = UpdateP(); // 誤差共分散行列を更新
#endif // TIME_UPDATE_DEBUG

    for (int i = 0; i < 2; i++)
    {
      pre_visible_gnss_nums_.at(i) = visible_gnss_nums_.at(i);
    }
    visible_gnss_nums_.clear();
    gnss_observations_.clear();

    PBD_GnssObservation main_observation = main_observation_;
    PBD_GnssObservation target_observation = target_observation_;
    gnss_observations_.push_back(main_observation); // 参照しない方がいい．
    gnss_observations_.push_back(target_observation);

    visible_gnss_nums_.push_back(main_observation_.GetNowVisibleGnssNum());
    visible_gnss_nums_.push_back(target_observation_.GetNowVisibleGnssNum());
    int n_main = pre_visible_gnss_nums_.at(0);
    Eigen::MatrixXd P_main = P_.topLeftCorner(NUM_SINGLE_STATE + n_main, NUM_SINGLE_STATE + n_main);
    Eigen::MatrixXd Q_main = Q_.topLeftCorner(NUM_SINGLE_STATE + n_main, NUM_SINGLE_STATE + n_main);
    int n_target = pre_visible_gnss_nums_.at(1);
    Eigen::MatrixXd P_target = P_.bottomRightCorner(NUM_SINGLE_STATE + n_target, NUM_SINGLE_STATE + n_target);
    Eigen::MatrixXd Q_target = Q_.bottomRightCorner(NUM_SINGLE_STATE + n_target, NUM_SINGLE_STATE + n_target);

    // 共通衛星見つける
    FindCommonObservedGnss(std::make_pair(0, 1));
    visible_gnss_nums_.push_back(common_observed_gnss_sat_id.size());

    UpdateBiasForm(0, x_est_main, P_main, Q_main);
    UpdateBiasForm(1, x_est_target, P_target, Q_target);
    const int new_size_all = P_main.rows() + P_target.rows();
    P_ = Eigen::MatrixXd::Zero(new_size_all, new_size_all);
    Q_ = Eigen::MatrixXd::Zero(new_size_all, new_size_all);
    int n_main_new = visible_gnss_nums_.at(0);
    P_.topLeftCorner(NUM_SINGLE_STATE + n_main_new, NUM_SINGLE_STATE + n_main_new) = P_main;
    Q_.topLeftCorner(NUM_SINGLE_STATE + n_main_new, NUM_SINGLE_STATE + n_main_new) = Q_main;
    int n_target_new = visible_gnss_nums_.at(1);
    P_.bottomRightCorner(NUM_SINGLE_STATE + n_target_new, NUM_SINGLE_STATE + n_target_new) = P_target;
    Q_.bottomRightCorner(NUM_SINGLE_STATE + n_target_new, NUM_SINGLE_STATE + n_target_new) = Q_target;

    KalmanFilter();

#ifndef TIME_UPDATE_DEBUG
    InitializePhi();
#endif // TIME_UPDATE_DEBUG

    // ここで観測情報を次用に更新する．
    main_observation_.UpdateInfoAfterObserved();
    target_observation_.UpdateInfoAfterObserved();
  }

  //log output
  if(abs(elapsed_time - tmp_log*log_step_time) < step_time/2.0){
    libra::Vector<3> sat_position_main = main_dynamics_.GetPosition_i();
    libra::Vector<3> sat_velocity_main = main_dynamics_.GetOrbit().GetSatVelocity_i();
    libra::Vector<3> sat_position_target = target_dynamics_.GetPosition_i();
    libra::Vector<3> sat_velocity_target = target_dynamics_.GetOrbit().GetSatVelocity_i();

    // ECIでの真値（位置，クロックバイアス，速度）を残す．
    for(int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << sat_position_main[i] << ","; // r_m_true
    ofs << std::fixed << std::setprecision(precision) << gnss_observations_.at(0).GetReceiver()->GetClockBias() << ","; // t_m_true
    for(int i = 0;i < 3;++i) ofs << std::fixed << std::setprecision(precision) << sat_velocity_main[i] << ","; // v_m_true

    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << sat_position_target[i] << ","; // r_t_true
    ofs << std::fixed << std::setprecision(precision) << gnss_observations_.at(1).GetReceiver()->GetClockBias() << ","; // t_t_true
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << sat_velocity_target[i] << ","; // v_t_true

    // 推定結果，ECIでの値を残す．
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_main.position(i) << ","; // r_m_est
    ofs << std::fixed << std::setprecision(precision) << x_est_main.clock(0) << ","; // t_m_est
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_main.velocity(i) << ","; // v_m_est
    for(int i = 0;i < 3;++i) ofs << std::fixed << std::setprecision(precision) << x_est_main.acceleration(i) << ","; // a_m_est

    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_target.position(i) << ","; // r_t_est
    ofs << std::fixed << std::setprecision(precision) << x_est_target.clock(0) << ","; //t_t_est
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_target.velocity(i) << ","; // v_t_est
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_target.acceleration(i) << ","; // a_t_est

    Eigen::Vector3d sat_pos_eci_main{ };
    Eigen::Vector3d sat_vel_eci_main{ };
    Eigen::Vector3d sat_pos_eci_target{ };
    Eigen::Vector3d sat_vel_eci_target{ };
    for (uint8_t i = 0; i < 3; i++)
    {
      sat_pos_eci_main(i) = sat_position_main[i];
      sat_pos_eci_target(i) = sat_position_target[i];
      sat_vel_eci_main(i) = sat_velocity_main[i];
      sat_vel_eci_target(i) = sat_velocity_target[i];
    }
    Eigen::Matrix3d trans_rtn_to_eci_main = TransRTN2ECI(sat_pos_eci_main, sat_vel_eci_main);
    Eigen::Matrix3d trans_eci_to_rtn_main = trans_rtn_to_eci_main.inverse();
    Eigen::Matrix3d trans_rtn_to_eci_target = TransRTN2ECI(sat_pos_eci_target, sat_vel_eci_target);
    Eigen::Matrix3d trans_eci_to_rtn_target = trans_rtn_to_eci_target.inverse();

    // RTNでの残差（位置，速度）を残す．
    Eigen::Vector3d res_pos_rtn_main{ };
    Eigen::Vector3d res_vel_rtn_main{ };
    Eigen::Vector3d res_pos_rtn_target{ };
    Eigen::Vector3d res_vel_rtn_target{ };
    res_pos_rtn_main = trans_eci_to_rtn_main * (x_est_main.position - sat_pos_eci_main);
    res_vel_rtn_main = trans_eci_to_rtn_main * (x_est_main.velocity - sat_vel_eci_main);
    res_pos_rtn_target = trans_eci_to_rtn_target * (x_est_target.position - sat_pos_eci_target);
    res_vel_rtn_target = trans_eci_to_rtn_target * (x_est_target.velocity - sat_vel_eci_target);
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << res_pos_rtn_main(i) << ","; // res_pos_m_rtn
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << res_vel_rtn_main(i) << ","; // res_vel_m_rtn
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << res_pos_rtn_target(i) << ","; // res_pos_t_rtn
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << res_vel_rtn_target(i) << ","; // res_vel_t_rtn

    for (int i = 0; i < NUM_GNSS_CH; ++i) ofs << std::fixed << std::setprecision(precision) << sat_info_.at(0).true_N(i) << ","; // N_true
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i < visible_gnss_nums_.at(0)) ofs << std::fixed << std::setprecision(precision) << x_est_main.ambiguity.N.at(i) << ","; // N_est
      else ofs << 0 << ",";
    }
    for (int i = 0; i < NUM_GNSS_CH; ++i) ofs << std::fixed << std::setprecision(precision) << sat_info_.at(1).true_N(i) << ","; // N_true
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i < visible_gnss_nums_.at(1)) ofs << std::fixed << std::setprecision(precision) << x_est_target.ambiguity.N.at(i) << ","; // N_est
      else ofs << 0 << ",";
    }

    const int non_visible_num_main = NUM_GNSS_CH - visible_gnss_nums_.at(0);
    Eigen::MatrixXd P_main = P_.topLeftCorner(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), NUM_SINGLE_STATE + visible_gnss_nums_.at(0));
    LogOutput_(ofs, P_main, NUM_SINGLE_STATE + visible_gnss_nums_.at(0), NUM_SINGLE_STATE_ALL);
    Eigen::MatrixXd P_target = P_.bottomRightCorner(NUM_SINGLE_STATE + visible_gnss_nums_.at(1), NUM_SINGLE_STATE + visible_gnss_nums_.at(1));
    LogOutput_(ofs, P_target, NUM_SINGLE_STATE + visible_gnss_nums_.at(1), NUM_SINGLE_STATE_ALL);

    // RTNでのcovariance(r, vのみ)
    TransECI2RTN_P(P_main, trans_eci_to_rtn_main);
    for (int i = 0; i < 3; i++) ofs << std::fixed << std::setprecision(precision) << P_main(i, i) << ","; // P_rtn_main (position)
    for (int i = 0; i < 3; i++) ofs << std::fixed << std::setprecision(precision) << P_main(4 + i, 4 + i) << ","; // P_rtn_main (velocity)
    TransECI2RTN_P(P_target, trans_eci_to_rtn_target);
    for (int i = 0; i < 3; i++) ofs << std::fixed << std::setprecision(precision) << P_target(i, i) << ","; // P_rtn_target (position)
    for (int i = 0; i < 3; i++) ofs << std::fixed << std::setprecision(precision) << P_target(4 + i, 4 + i) << ","; // P_rtn_target (velocity)

    // record visible gnss sat number
    // そもそもここでログをとるのが適切ではない．
    ofs << visible_gnss_nums_.at(0) << ",";
    ofs << visible_gnss_nums_.at(1) << ",";
    ofs << visible_gnss_nums_.at(2) << ",";
    // main observe gnss sat id
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i >= visible_gnss_nums_.at(0)) ofs << -1 << ",";
      else ofs << gnss_observations_.at(0).info_.now_observed_gnss_sat_id.at(i) << ",";
    }
    // target observe gnss sat id
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i >= visible_gnss_nums_.at(1)) ofs << -1 << ",";
      else ofs << gnss_observations_.at(1).info_.now_observed_gnss_sat_id.at(i) << ",";
    }

    // Q_
      Eigen::MatrixXd Q_main = Q_.topLeftCorner(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), NUM_SINGLE_STATE + visible_gnss_nums_.at(0));
    LogOutput_(ofs, Q_main, NUM_SINGLE_STATE + visible_gnss_nums_.at(0), NUM_SINGLE_STATE_ALL);
    Eigen::MatrixXd Q_target = Q_.bottomRightCorner(NUM_SINGLE_STATE + visible_gnss_nums_.at(1), NUM_SINGLE_STATE + visible_gnss_nums_.at(1));
    LogOutput_(ofs, Q_target, NUM_SINGLE_STATE + visible_gnss_nums_.at(1), NUM_SINGLE_STATE_ALL);

    // R_
    Eigen::MatrixXd R_gr_m = R_.topLeftCorner(visible_gnss_nums_.at(0), visible_gnss_nums_.at(0));
    LogOutput_(ofs, R_gr_m, visible_gnss_nums_.at(0), NUM_GNSS_CH);
    Eigen::MatrixXd R_gr_t = R_.block(visible_gnss_nums_.at(0), visible_gnss_nums_.at(0), visible_gnss_nums_.at(1), visible_gnss_nums_.at(1));
    LogOutput_(ofs, R_gr_t, visible_gnss_nums_.at(1), NUM_GNSS_CH);
    Eigen::MatrixXd R_sdcp = R_.bottomRightCorner(visible_gnss_nums_.at(2), visible_gnss_nums_.at(2));
    LogOutput_(ofs, R_sdcp, visible_gnss_nums_.at(2), NUM_GNSS_CH);

    // acc eci
    Eigen::Vector3d acc_m_i = trans_rtn_to_eci_main* x_est_main.acceleration; // [nm/s2]
    Eigen::Vector3d acc_t_i = trans_rtn_to_eci_target * x_est_target.acceleration; // [nm/s2]
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << acc_m_i(i) << ","; // a_m_i
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << acc_t_i(i) << ","; // a_t_i
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_main.acc_dist(i) << ","; // a_disturbance_m
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_target.acc_dist(i) << ","; // a_disturbance_t
    // acc rtn
    Eigen::Vector3d acc_dist_m_rtn = trans_eci_to_rtn_main* x_est_main.acc_dist; // [m/s2]
    Eigen::Vector3d acc_dist_t_rtn = trans_eci_to_rtn_target * x_est_target.acc_dist; // [m/s2]
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << acc_dist_m_rtn(i) << ","; // a_disturbance_m_rtn
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << acc_dist_t_rtn(i) << ","; // a_disturbance_t_rtn

    // azimuth elevation
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i < visible_gnss_nums_.at(0)) ofs << std::fixed << std::setprecision(precision) << gnss_observations_.at(0).GetGnssAzimuthDeg(i) << ","; // azimuth main
      else ofs << 0 << ",";
    }
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i < visible_gnss_nums_.at(0)) ofs << std::fixed << std::setprecision(precision) << gnss_observations_.at(0).GetGnssElevationDeg(i) << ","; // elevation main
      else ofs << 0 << ",";
    }
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i < visible_gnss_nums_.at(1)) ofs << std::fixed << std::setprecision(precision) << gnss_observations_.at(1).GetGnssAzimuthDeg(i) << ","; // azimuth target
      else ofs << 0 << ",";
    }
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i < visible_gnss_nums_.at(1)) ofs << std::fixed << std::setprecision(precision) << gnss_observations_.at(1).GetGnssElevationDeg(i) << ","; // elevation target
      else ofs << 0 << ",";
    }

    ofs << std::endl;
  }

  return;
}

// もう少し汎用性の高い形にする．
static void LogOutput_(std::ofstream& ofs, const Eigen::MatrixXd& M, const int size, const int max_size)
{
  for (int i = 0; i < max_size; ++i)
  {
    if (i < size) ofs << std::fixed << std::setprecision(precision) << M(i, i) << ",";
    else ofs << 0 << ",";
  }
}

void PBD_dgps::InitializePhi(void)
{
  Phi_ = Eigen::Matrix<double, NUM_STATE, NUM_STATE>::Identity();
  Phi_(3, 3) = 0.0;
  Phi_(NUM_SINGLE_STATE + 3, NUM_SINGLE_STATE + 3) = 0.0;
}

void PBD_dgps::TransECI2RTN_P(Eigen::MatrixXd& P, Eigen::Matrix3d trans_eci_to_rtn)
{
  Eigen::Matrix3d P_pos = P.topLeftCorner(3, 3);
  Eigen::Matrix3d P_vel = P.block(4, 4, 3, 3);

  P.topLeftCorner(3, 3) = trans_eci_to_rtn * P_pos * trans_eci_to_rtn.transpose();
  P.block(4, 4, 3, 3) = trans_eci_to_rtn * P_vel * trans_eci_to_rtn.transpose();
}


void PBD_dgps::InitLogTable(void)
{

}

void PBD_dgps::OrbitPropagation()
{
  // Phiは観測更新間の状態遷移行列．
  Eigen::MatrixXd Phi_main = Phi_.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE);
  Eigen::MatrixXd Phi_target = Phi_.bottomRightCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE);

  RK4(x_est_main.position, x_est_main.velocity, x_est_main.acceleration, x_est_main.acc_dist, Phi_main);
  Phi_.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Phi_main;

  RK4(x_est_target.position, x_est_target.velocity, x_est_target.acceleration, x_est_target.acc_dist, Phi_target);
  Phi_.bottomRightCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Phi_target;

#ifdef REDUCED_DYNAMIC
  // acceleration ここでは q = sigma_acc_process
  Eigen::MatrixXd Phi_a = CalculatePhi_a(observe_step_time); // ここもチューニングできないと意味ないのでは？
  double phi = Phi_a(0, 0);
  // double phi_t = Phi_a(1, 1);
  // double phi_n = Phi_a(2, 2);
  // x_est_main.acceleration = Phi_a * x_est_main.acceleration;
  // x_est_target.acceleration = Phi_a * x_est_target.acceleration;

#endif // REDUCED_DYNAMIC

  // cdt
#ifdef CLOCK_IS_RANDOM_WALK
  std::normal_distribution<> cdt_process_noise(0.0, sigma_cdt_process*sqrt(step_time/tau_cdt));
  std::normal_distribution<> cdt_process_noise(0.0, sigma_cdt_process);
  x_est_main.clock(0) = cdt_process_noise(mt);
  x_est_target.clock(0) = cdt_process_noise(mt);
#endif // CLOCK_IS_RANDOM_WALK

#ifdef TIME_UPDATE_DEBUG
  // Pの伝搬はここでしないとモデル化誤差がでかいときに死んでしまうのかもしれない．
  P_ = UpdateP();
  InitializePhi();
#endif // TIME_UPDATE_DEBUG
}

void PBD_dgps::RK4(Eigen::Vector3d& position, Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration, Eigen::Vector3d& acc_dist, Eigen::MatrixXd& Phi)
{
  Eigen::Vector3d k0 = PositionDifferential(velocity);
  Eigen::Vector3d l0 = VelocityDifferential(position, velocity, acceleration, acc_dist);
  Eigen::Vector3d m0 = -acceleration / tau_a;
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> n0 = CalculateJacobian(position, velocity); //*Phi;

  Eigen::Vector3d tmp_position = position + k0 * step_time / 2.0;
  Eigen::Vector3d tmp_velocity = velocity + l0 * step_time / 2.0;
  Eigen::Vector3d tmp_acceleration = acceleration + m0 * step_time / 2.0;
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> tmp_Phi = Phi + n0 * step_time / 2.0;
  Eigen::Vector3d k1 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l1 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Vector3d m1 = -tmp_acceleration / tau_a;
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> n1 = CalculateJacobian(tmp_position, tmp_velocity); //*tmp_Phi;

  tmp_position = position + k1 * step_time / 2.0;
  tmp_velocity = velocity + l1 * step_time / 2.0;
  tmp_acceleration = acceleration + m1 * step_time / 2.0;
  tmp_Phi = Phi + n1 * step_time / 2.0;
  Eigen::Vector3d k2 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l2 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Vector3d m2 = -tmp_acceleration / tau_a;
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> n2 = CalculateJacobian(tmp_position, tmp_velocity); // *tmp_Phi;

  tmp_position = position + k2 * step_time;
  tmp_velocity = velocity + l2 * step_time;
  tmp_acceleration = acceleration + m2 * step_time;
  tmp_Phi = Phi + n2 * step_time;
  Eigen::Vector3d k3 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l3 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Vector3d m3 = -tmp_acceleration / tau_a;
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> n3 = CalculateJacobian(tmp_position, tmp_velocity); // *tmp_Phi;

  position     += step_time * (k0 + 2.0 * k1 + 2.0 * k2 + k3) / 6.0;
  velocity     += step_time * (l0 + 2.0 * l1 + 2.0 * l2 + l3) / 6.0;
  acceleration += step_time * (m0 + 2.0 * m1 + 2.0 * m2 + m3) / 6.0; // 加速度に関しては解析モデルを使った方がいい気がする．
  Phi          += step_time * (n0 + 2.0 * n1 + 2.0 * n2 + n3) / 6.0;
}

Eigen::Vector3d PBD_dgps::PositionDifferential(const Eigen::Vector3d& velocity) const
{
  return velocity;
}

Eigen::Vector3d PBD_dgps::VelocityDifferential(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration, Eigen::Vector3d& acc_dist) const
{
  double r = position.norm();
  double v = velocity.norm();

  Eigen::Vector3d acc_2body = - mu_e * position / pow(r, 3.0); // 2体の重力項

  // libra::Vector<3> position_eci;
  // for (uint8_t i = 0; i < 3; i++) position_eci[i] = position(i);
  // libra::Vector<3> position_ecef = trans_eci_to_ecef_ * position_eci;

  // double x = position_ecef[0];
  // double y = position_ecef[1];
  // double z = position_ecef[2];
  // double tmp_J2_coefficient = 3.0/2.0*mu_e*J2_const*pow(Earth_Radius, 2.0)/pow(r, 4.0); // J2項の係数
  // libra::Vector<3> acc_j2_ecef;
  // acc_j2_ecef[0] = -(tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0))) * (x / r);
  // acc_j2_ecef[1] = -(tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0))) * (y / r);
  // acc_j2_ecef[2] = -(tmp_J2_coefficient*(3.0 - 5.0*pow(z/r, 2.0))) * (z / r);

  // libra::Matrix<3, 3> trans_ecef_to_eci =  libra::transpose(trans_eci_to_ecef_);
  // libra::Vector<3> acc_j2_eci = trans_ecef_to_eci * acc_j2_ecef;
  // for (uint8_t i = 0; i < 3; i++) acc_dist(i) = acc_j2_eci[i];

#ifndef GEOP_DEBUG
  AddGeoPotentialDisturbance(position, acc_dist);
#endif // GEOP_DEBUG

#ifdef AIR_DRAG_ON
  acc_dist -= Cd*v*velocity; // -Cd*V^2*(Vi/V) 大気抵抗
#endif // AIR_DRAG_ON

  Eigen::MatrixXd empirical_acc(3, 1);
  empirical_acc.block(0, 0, 3, 1) = acceleration; // RTN
  Eigen::Vector3d empirical_acc_eci = TransRTN2ECI(position, velocity) * empirical_acc;

  Eigen::Vector3d all_acceleration = acc_2body + acc_dist +  empirical_acc_eci*PBD_DGPS_kConvNm2m;
  return all_acceleration; // m/s2
}

void PBD_dgps::AddGeoPotentialDisturbance(const Eigen::Vector3d& position, Eigen::Vector3d& acc_dist) const
{
  libra::Vector<3> position_eci;
  for (uint8_t i = 0; i < 3; i ++) position_eci[i] = position(i);

  libra::Vector<3> acc_geop;
  acc_geop = geo_potential_->CalcAccelerationECI(position_eci, trans_eci_to_ecef_);
  for (uint8_t i = 0; i < 3; i ++) acc_dist(i) = acc_geop[i];
}

Eigen::MatrixXd PBD_dgps::UpdateP(void)
{
  int num_state_main = NUM_SINGLE_STATE;
  int num_state_target = NUM_SINGLE_STATE;
  int n_main = 0;
  int n_target = 0;
  if (visible_gnss_nums_.at(0) != 0)
  {
    n_main = visible_gnss_nums_.at(0);
    num_state_main += n_main;
    n_target = visible_gnss_nums_.at(1);
    num_state_target += n_target;
  }
  const int num_state_all = num_state_main + num_state_target;
  Eigen::MatrixXd Phi_all = Eigen::MatrixXd::Identity(num_state_all, num_state_all);

  Phi_all.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Phi_.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE);
  Phi_all.block(num_state_main, num_state_main, NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Phi_.bottomRightCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE);

  // clock
#ifdef CLOCK_IS_RANDOM_WALK
  Phi_all(3, 3) = 1.0; // random walkなのでΦは1
  Phi_all(num_state_main + 3, num_state_main + 3) = 1.0;
#endif // CLOCK_IS_RANDOM_WALK
#ifdef CLOCK_IS_WHITE_NOISE
  // ここは普通に伝搬するのが正解そう?
  Phi_all(3, 3) = 0.0;
  Phi_all(num_state_main + 3, num_state_main + 3) = 0.0;
#endif // CLOCK_IS_WHITE_NOISE

  // a
#ifdef REDUCED_DYNAMIC
  Phi_all.block(7, 7, 3, 3) = CalculatePhi_a(observe_step_time);
  Phi_all.block(num_state_main + 7, num_state_main + 7, 3, 3) = CalculatePhi_a(observe_step_time);
#endif // REDUCED_DYNAMIC

#ifdef AKF
  // ここでは更新しない．
#else
  InitializeQ();
#endif // AKF

  Eigen::MatrixXd res = Phi_all * P_ * Phi_all.transpose() + Q_; // ここを毎回加えるのはダメそう．

  return res;
}

Eigen::MatrixXd PBD_dgps::CalculateJacobian(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const
{
  double r = position.norm(); // [m]
  double v = velocity.norm(); // [m/s]

  libra::Vector<3> position_eci;
  for (uint8_t i = 0; i < 3; i++) position_eci[i] = position(i);
  libra::Vector<3> position_ecef = trans_eci_to_ecef_ * position_eci;
  double x_i = position_eci[0]; double y_i = position_eci[1]; double z_i = position_eci[2];
  double x_ef = position_ecef[0]; double y_ef = position_ecef[1]; double z_ef = position_ecef[2];
  double vx = velocity(0); double vy = velocity(1); double vz = velocity(2);


  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(NUM_SINGLE_STATE, NUM_SINGLE_STATE);
  // (r, v)
  A(0, 4) = 1.0; A(1, 5) = 1.0; A(2, 6) = 1.0;
  // (v, r)
  double r3 = pow(r, 3);
  double r5 = pow(r, 5);
  A(4, 0) = mu_e / r3 * (3.0 * pow(x_i / r, 2.0) - 1); A(4, 1) = mu_e / r5 * 3.0 * x_i*y_i; A(4, 2) = mu_e / r5 * 3.0 * x_i*z_i;
  A(5, 0) = mu_e * 3.0 * y_i*x_i / r5; A(5, 1) = mu_e / r3 * (3.0 * pow(y_i / r, 2.0) - 1); A(5, 2) = mu_e * 3.0 * y_i*z_i / r5;
  A(6, 0) = mu_e * 3.0 * z_i*x_i / r5; A(6, 1) = mu_e * 3.0 * z_i*y_i / r5; A(6, 2) = mu_e / r3 * (3.0 * pow(z_i / r, 2.0) - 1);

  #ifndef GEOP_DEBUG
  // J2
  double J2_coefficient = 3.0 / 2.0 * mu_e * J2_const * pow(Earth_Radius, 2.0);
  libra::Matrix<3, 3> A_j2_ecef(0);
  A_j2_ecef[0][0] = - J2_coefficient / pow(r, 5.0) * (1.0 - 5.0 * (pow(x_ef, 2.0) + pow(z_ef, 2.0)) / pow(r, 2.0) + 35.0 * pow(x_ef * z_ef, 2.0) / pow(r, 4.0));
  A_j2_ecef[0][1] = - J2_coefficient * x_ef * y_ef / pow(r, 7.0) * (-5.0 + 35.0 * pow(z_ef / r, 2.0));
  A_j2_ecef[0][2] = - J2_coefficient * x_ef * z_ef / pow(r, 7.0) * (-15.0 + 35.0 * pow(z_ef / r, 2.0));
  A_j2_ecef[1][0] = - J2_coefficient * x_ef * y_ef / pow(r, 7.0) * (-5.0 + 35.0 * pow(z_ef / r, 2.0));
  A_j2_ecef[1][1] = - J2_coefficient / pow(r, 5.0) * (1.0 - 5.0 * (pow(y_ef, 2.0) + pow(z_ef, 2.0)) / pow(r, 2.0) + 35.0 * pow(y_ef * z_ef, 2.0) / pow(r, 4.0));
  A_j2_ecef[1][2] = - J2_coefficient * y_ef * z_ef / pow(r, 7.0) * (-15.0 + 35.0 * pow(z_ef / r, 2.0));
  A_j2_ecef[2][0] = - J2_coefficient * x_ef * z_ef / pow(r, 7.0) * (-15.0 + 35.0 * pow(z_ef / r, 2.0));
  A_j2_ecef[2][1] = - J2_coefficient * y_ef * z_ef / pow(r, 7.0) * (-15.0 + 35.0 * pow(z_ef / r, 2.0));
  A_j2_ecef[2][2] = - J2_coefficient / pow(r, 5.0) * (3.0 - 30.0 * pow(z_ef/r, 2.0) + 35.0 * pow(z_ef/r, 4.0));

  libra::Matrix<3, 3> A_j2_eci = libra::transpose(trans_eci_to_ecef_) * A_j2_ecef * trans_eci_to_ecef_;

  for (uint8_t i = 0; i < 3; i++)
  {
    for (uint8_t j = 0; j < 3; j++) A(i + 4, j) += A_j2_eci[i][j];
  }
  #endif // GEOP_DEBUG

#ifdef AIR_DRAG_ON
  // 空気抵抗分
  A(4, 4) = -Cd * (vx * vx / v + v);    A(4, 5) = -Cd * vx * vy / v;    A(4, 6) = -Cd * vx * vz / v;
  A(5, 4) = -Cd * vx * vy / v;    A(5, 5) = -Cd * (vy * vy / v + v);    A(5, 6) = -Cd * vy * vz / v;
  A(6, 4) = -Cd * vx * vz / v;    A(6, 5) = -Cd * vy * vz / v;    A(6, 6) = -Cd * (vz * vz / v + v);
#endif // AIR_DRAG_ON

  // (v, a)
#ifdef REDUCED_DYNAMIC
  Eigen::Matrix3d rtn2eci = TransRTN2ECI(position, velocity);
  A.block(4, 7, 3, 3) = PBD_DGPS_kConvNm2m*rtn2eci;
#endif // REDUCED_DYNAMIC

  // これはもしCdも推定しているならいる．
  //A(4,10) = -v*vx;    A(5,10) = -v*vy;    A(6,10) = -v*vz;
  return A;
};

// これはlibra::Vectorにした方がいいかもしれん．
Eigen::Matrix3d PBD_dgps::TransRTN2ECI(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const
{
  Eigen::Vector3d r_eci = position.normalized();
  Eigen::Vector3d n_eci = position.cross(velocity);
  n_eci.normalize();
  Eigen::Vector3d t_eci = n_eci.cross(r_eci);
  t_eci.normalize();

  Eigen::MatrixXd RTN2ECI(3,3);
  RTN2ECI.block(0, 0, 3, 1) = r_eci;
  RTN2ECI.block(0, 1, 3, 1) = t_eci;
  RTN2ECI.block(0, 2, 3, 1) = n_eci;
  return RTN2ECI;
};

void PBD_dgps::InitializeQ(void)
{
  const int n_main = visible_gnss_nums_.at(0);
  const int n_target = visible_gnss_nums_.at(1);
  const int num_state_all = NUM_STATE + n_main + n_target;
  const int num_state_main = NUM_SINGLE_STATE + n_main;
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_state_all, 8); // aとclock
  B(3, 0) = 1.0;
  B(num_state_main + 3, 4) = 1.0;
#ifdef REDUCED_DYNAMIC
  B.block(7, 1, 3, 3) = Eigen::Matrix3d::Identity();
  B.block(num_state_main + 7, 5, 3, 3) = Eigen::Matrix3d::Identity();
#endif // REDUCED_DYNAMIC

  Eigen::MatrixXd Q_at = CalculateQ_at();
  // 観測するたびにNの部分を初期化？

  // Q_ = BQ_atB^t
  Q_ = B * Q_at * B.transpose();

  // add process noise for r and v
  Q_.block(0, 0, 3, 3) = pow(sigma_r_process, 2.0) * Eigen::Matrix3d::Identity(); // * pow(step_time, 2.0);
  Q_.block(4, 4, 3, 3) = pow(sigma_v_process, 2.0) * Eigen::Matrix3d::Identity(); // * pow(step_time, 2.0);
  Q_.block(num_state_main, num_state_main, 3, 3) = pow(sigma_r_process, 2.0) * Eigen::Matrix3d::Identity(); // * pow(step_time, 2.0);
  Q_.block(num_state_main + 4, num_state_main + 4, 3, 3) = pow(sigma_v_process, 2.0) * Eigen::Matrix3d::Identity(); // * pow(step_time, 2.0);
#ifdef TIME_UPDATE_DEBUG
  Q_.block(0, 0, 3, 3) *= step_time / observe_step_time;
  Q_.block(4, 4, 3, 3) *= step_time / observe_step_time;
  Q_.block(num_state_main, num_state_main, 3, 3) *= step_time / observe_step_time;
  Q_.block(num_state_main + 4, num_state_main + 4, 3, 3) *= step_time / observe_step_time;
#endif // TIME_UPDATE_DEBUG

#ifndef N_DEBUG
  // N process
  Q_.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, n_main, n_main) = pow(sigma_N_process, 2.0) * Eigen::MatrixXd::Identity(n_main, n_main); // * pow(step_time, 2.0);
  Q_.block(NUM_SINGLE_STATE + num_state_main, NUM_SINGLE_STATE + num_state_main, n_target, n_target) = pow(sigma_N_process, 2.0) * Eigen::MatrixXd::Identity(n_target, n_target); // * pow(step_time, 2.0);
#ifdef TIME_UPDATE_DEBUG
  Q_.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, n_main, n_main) *= step_time / observe_step_time;
  Q_.block(NUM_SINGLE_STATE + num_state_main, NUM_SINGLE_STATE + num_state_main, n_target, n_target) *= step_time / observe_step_time;
#endif // TIME_UPDATE_DEBUG
#endif // N_DEBUG
}

// Process noiseのvarianceを計算．Bを使う形に修正．
Eigen::MatrixXd PBD_dgps::CalculateQ_at(void)
{
  Eigen::MatrixXd Q_at = Eigen::MatrixXd::Zero(8, 8);
  double q_acc_r = 0.0f;
  double q_acc_t = 0.0f;
  double q_acc_n = 0.0f;

#ifdef REDUCED_DYNAMIC
  double phi = exp(-observe_step_time / tau_a); // ここはobserved_timeとどっちなのか？
  q_acc_r = pow(sigma_acc_r_process, 2.0) * (1 - pow(phi, 2.0));
  q_acc_t = pow(sigma_acc_t_process, 2.0) * (1 - pow(phi, 2.0));
  q_acc_n = pow(sigma_acc_n_process, 2.0) * (1 - pow(phi, 2.0));
#endif // REDUCED_DYNAMIC

  double q_cdt;
#ifdef CLOCK_IS_WHITE_NOISE
  q_cdt = pow(sigma_cdt_process, 2.0); //  * step_time
#endif // CLOCK_IS_WHITE_NOISE
#ifdef CLOCK_IS_RANDOM_WALK
  q_cdt = pow(sigma_cdt_process, 2.0) * (observe_step_time / tau_cdt);
#endif // CLOCK_IS_RANDOM_WALK

  Q_at(0, 0) = q_cdt;
  Q_at(1, 1) = q_acc_r; Q_at(2, 2) = q_acc_t; Q_at(3, 3) = q_acc_n;

  Q_at(4, 4) = q_cdt;
  Q_at(5, 5) = q_acc_r; Q_at(6, 6) = q_acc_t; Q_at(7, 7) = q_acc_n;

  return Q_at;
}

Eigen::MatrixXd PBD_dgps::CalculatePhi_a(const double dt)
{
  double phi = exp(-dt / tau_a); // constなので計算毎回する意味ない．
  Eigen::MatrixXd Phi_a = Eigen::MatrixXd::Identity(3, 3);
#ifdef REDUCED_DYNAMIC
  Phi_a *= phi;
#endif // REDUCED_DYNAMIC

  return Phi_a;
};

void PBD_dgps::KalmanFilter()
{
  // gnss_observed_hoge の観測される時間はどうなっているのか？
  int n_main = visible_gnss_nums_.at(0);
  int n_target = visible_gnss_nums_.at(1);
  int n_common = visible_gnss_nums_.at(2);
  const int num_observables = n_main + n_target + n_common;
  const int num_main_state_all = NUM_SINGLE_STATE + n_main;
  const int num_state_all = NUM_STATE + n_main + n_target;

  // GRAPHIC*2 + SDCPにする．
  Eigen::VectorXd z = Eigen::VectorXd::Zero(num_observables); //観測ベクトル
  Eigen::VectorXd h_x = Eigen::VectorXd::Zero(num_observables); // 観測モデル行列

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(num_observables, num_state_all); //観測行列（dhx/dx）
  Eigen::VectorXd R_V = Eigen::MatrixXd::Zero(num_observables, 1); // 観測誤差共分散．
  UpdateObservations(z, h_x, H, R_V);

// #ifdef AKF
//   // Rのresizeが必要．
// #else
  R_ = R_V.asDiagonal();
// #endif

  Eigen::MatrixXd hmh = H * P_ * H.transpose();

  Eigen::MatrixXd tmp = R_ + hmh; // (observation_num, observation_num)
  // 観測量がないところは0にし直すみたいな加工が必要かもしれない．数字の桁数で打ち切りみたいな形にできればいい．
  //カルマンゲイン
  Eigen::MatrixXd K = CalculateK(H, tmp);

  Eigen::VectorXd x_predict = Eigen::VectorXd(num_state_all);
  x_predict.topRows(3) = x_est_main.position;
  x_predict.block(3, 0, 1, 1) = x_est_main.clock;
  x_predict.block(4, 0, 3, 1)  = x_est_main.velocity;
#ifdef REDUCED_DYNAMIC
  x_predict.block(7, 0, 3, 1) = x_est_main.acceleration;
  x_predict.block(num_main_state_all + 7, 0, 3, 1) = x_est_target.acceleration;
#endif // REDUCED_DYNAMIC
  x_predict.block(NUM_SINGLE_STATE, 0, n_main, 1) = Eigen::Map<Eigen::VectorXd>(&x_est_main.ambiguity.N[0], n_main);
  //x_predict(10) = Cd;
  x_predict.block(num_main_state_all, 0, 3, 1) = x_est_target.position;
  x_predict.block(num_main_state_all + 3, 0, 1, 1) = x_est_target.clock;
  x_predict.block(num_main_state_all + 4, 0, 3, 1) = x_est_target.velocity;
  x_predict.bottomRows(n_target) = Eigen::Map<Eigen::VectorXd>(&x_est_target.ambiguity.N[0], n_target);

  // innovationの記号を何にするかは要検討
  Eigen::VectorXd E_pre = z - h_x;
  // まずアンテナ位置に変換
  Eigen::VectorXd x_ant_predict = x_predict;
  x_ant_predict.topRows(3) = ConvCenterOfMassToReceivePos(x_predict.topRows(3), antenna_pos_b_.at(0), main_dynamics_);
  x_ant_predict.block(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), 0, 3, 1) = ConvCenterOfMassToReceivePos(x_predict.block(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), 0, 3, 1), antenna_pos_b_.at(1), target_dynamics_);
  Eigen::VectorXd x_update = x_ant_predict + K * E_pre;
  // Eigen::VectorXd x_update = x_predict + K * E_pre;
  // 重心位置に戻す．
  x_update.topRows(3) = ConvReceivePosToCenterOfMass(x_update.topRows(3), antenna_pos_b_.at(0), main_dynamics_);
  x_update.block(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), 0, 3, 1) = ConvReceivePosToCenterOfMass(x_update.block(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), 0, 3, 1), antenna_pos_b_.at(1), target_dynamics_);

 //更新
  x_est_main.position = x_update.topRows(3);
  x_est_main.clock = x_update.block(3, 0, 1, 1);
  x_est_main.velocity = x_update.block(4, 0, 3, 1);

  x_est_target.position = x_update.block(num_main_state_all, 0, 3, 1);
  x_est_target.clock = x_update.block(num_main_state_all + 3, 0, 1, 1);
  x_est_target.velocity = x_update.block(num_main_state_all + 4, 0, 3, 1);
#ifdef REDUCED_DYNAMIC
  x_est_main.acceleration = x_update.block(7, 0, 3, 1);
  x_est_target.acceleration = x_update.block(num_main_state_all + 7, 0, 3, 1);
#endif // REDUCED_DYNAMIC
#ifndef N_DEBUG
  for (int i = 0; i < NUM_GNSS_CH; ++i)
  {
    if (i < n_main) x_est_main.ambiguity.N.at(i) = x_update(NUM_SINGLE_STATE + i);
    if (i < n_target) x_est_target.ambiguity.N.at(i) = x_update(num_main_state_all + NUM_SINGLE_STATE + i);
  }
#else
    // 一旦整数不定性は0とする．
    x_est_main.ambiguity.N = std::vector<double>(NUM_GNSS_CH, 0);
    x_est_target.ambiguity.N = std::vector<double>(NUM_GNSS_CH, 0);
#endif // N_DEBUG

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(num_state_all, num_state_all);
  tmp = (I - K*H);
  P_ = tmp*P_*tmp.transpose() + K*R_*K.transpose();
  // P_ = tmp*P_;

#ifdef AKF
  // residual
  // GRAPHIC*2 + SDCPにする．
  Eigen::VectorXd h_x_post = Eigen::VectorXd::Zero(num_observables); // 観測モデル行列

  // clear
  for(int i = 0; i < 2; ++i) ClearGnssObserveModels(gnss_observed_models_.at(i));
  UpdateObservations(z, h_x_post, H, R_V);

  Eigen::VectorXd E_post = z - h_x_post;
  R_ = alpha * R_ + (1 - alpha) * (E_post*E_post.transpose() + H*P_*H.transpose()); // residual based R_-adaptation

  Eigen::MatrixXd Q_dash = alpha * Q_ + (1 - alpha) * K * E_pre * (K * E_pre).transpose(); // Innovation based Q_-adaptation

  DynamicNoiseScaling(Q_dash, H);
  // これしてるから，絶対軌道精度の影響（つまり残差の大きさ）ではなくて，収束してしまう？
/*
  Q_ = Eigen::MatrixXd::Zero(num_state_all, num_state_all);
  // a, cdtだけもらう．
  Q_(3, 3) = Q_dash(3, 3);
  Q_.block(7, 7, 3, 3) = Q_dash.block(7, 7, 3, 3);
  Q_(num_main_state_all + 3, num_main_state_all + 3) = Q_dash(num_main_state_all + 3, num_main_state_all + 3);
  Q_.block(num_main_state_all + 7, num_main_state_all + 7, 3, 3) = Q_dash.block(num_main_state_all + 7, num_main_state_all + 7, 3, 3);
*/
#endif // AKF

  // TODO: Nについて，収束している？ものは0に落とす．

  // これが必要になる時点で何かバグってそう．
  // for (int i = 0; i < NUM_GNSS_CH; ++i)
  // {
  //   int N_offset;
  //   if (i < n_main)
  //   {
  //     N_offset = NUM_SINGLE_STATE + i;
  //     if (P_(N_offset, N_offset) < 0) P_(N_offset, N_offset) = pow(sigma_N_ini, 2.0); // sqrt(pow(P_(N_offset, N_offset), 2.0));
  //   }
  //   if (i < n_target)
  //   {
  //     N_offset = NUM_STATE + n_main + i;
  //     if (P_(N_offset, N_offset) < 0) P_(N_offset, N_offset) = pow(sigma_N_ini, 2.0); // sqrt(pow(P_(N_offset, N_offset), 2.0));
  //   }
  // }

  // others
  // for (int i = 0; i < NUM_SINGLE_STATE; ++i)
  // {
  //   int offset = i;
  //   if (P_(offset, offset) < 1e-6) P_(offset, offset) = sqrt(pow(P_(offset, offset), 2.0));
  //   offset = NUM_SINGLE_STATE_ALL + i;
  //   if (P_(offset, offset) < 1e-6) P_(offset, offset) = sqrt(pow(P_(offset, offset), 2.0));
  // }

  // TODO: make the Double Differential ambiguity

  // IAR

#ifdef LAMBDA_DEBUG
  Eigen::MatrixXd Q_a_main = P_.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, NUM_GNSS_CH, NUM_GNSS_CH);
  Eigen::VectorXd a_main = x_est_main.N; // ambiguity
  // non-zeroに限定する．
  int num_of_obsetved_gnss =  gnss_observations_.at(0).GetVisibleGnssNum();
  PBD_Lambda lambda(Q_a_main.topLeftCorner(num_of_obsetved_gnss, num_of_obsetved_gnss), a_main.topRows(num_of_obsetved_gnss), 2);
  lambda.Solve();
#endif // LAMBDA_DEBUG

  if (!std::isfinite(x_est_main.position(0)))
  {
    std::cout << "inf or nan" << std::endl;
    abort();
  }

  // clear
  for(int i = 0; i < 2; ++i)
  {
    ClearGnssObserveModels(gnss_observed_models_.at(i));
  }
  return;
}

void PBD_dgps::ClearGnssObserveModels(GnssObserveModel& observed_model)
{
  observed_model.geometric_range.clear();
  observed_model.pseudo_range_model.clear();
  observed_model.carrier_phase_model.clear();
}

Eigen::MatrixXd PBD_dgps::CalculateK(Eigen::MatrixXd H, Eigen::MatrixXd S)
{
  const int n_main = visible_gnss_nums_.at(0);
  const int n_target = visible_gnss_nums_.at(1);
  const int n_common = visible_gnss_nums_.at(2);
  Eigen::MatrixXd PHt = P_ * H.transpose();
  // ResizeS(S, n_main, n_target, n_common); // 固定サイズの時に必要だったもの
  // ResizeMHt(MHt, n_main, n_target, n_common);
#if CHOLESKY
  Eigen::MatrixXd ST = S.transpose();
  Eigen::LDLT<Eigen::MatrixXd> LDLTOftmpT(ST);
  Eigen::MatrixXd KT = LDLTOftmpT.solve(PHt.transpose());
  Eigen::MatrixXd K = KT.transpose();
#elif QR
  Eigen::MatrixXd ST = S.transpose();
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QROftmpT(ST);
  Eigen::MatrixXd KT = QROftmpT.solve(PHt.transpose());
  Eigen::MatrixXd K = KT.transpose();
#else
  Eigen::MatrixXd K = PHt * S.inverse();
#endif

  return K;
};

// TODO: これをVectorとか他のファイルでも使えるように修正する．
void PBD_dgps::RemoveRows(Eigen::MatrixXd& matrix, unsigned int begin_row, unsigned int end_row)
{
  // FIXME: 全消しになった時にバグる．0行x列という行列は作れないので．
  // 消していくとずれるから後ろから
  for (int row = end_row; row >= begin_row; --row) {
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (row < numRows)
      matrix.block(row, 0, numRows - row, numCols) = matrix.bottomRows(numRows - row);

    matrix.conservativeResize(numRows, numCols);
  }
}

void PBD_dgps::RemoveColumns(Eigen::MatrixXd& matrix, unsigned int begin_col, unsigned int end_col)
{
  // 消していくとずれるから後ろから
  for (int col = end_col; col >= begin_col; --col) {

    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols() - 1;

  if (col < numCols)
    matrix.block(0, col, numRows, numCols - col) = matrix.rightCols(numCols - col);

  matrix.conservativeResize(numRows, numCols);
  }
};

void PBD_dgps::ResizeS(Eigen::MatrixXd& S, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c)
{
  // ここも後ろから．
  RemoveRows(S, 2*NUM_GNSS_CH + observe_gnss_c, 3*NUM_GNSS_CH - 1);
  RemoveColumns(S, 2*NUM_GNSS_CH + observe_gnss_c, 3*NUM_GNSS_CH - 1);
  RemoveRows(S, NUM_GNSS_CH + observe_gnss_t, 2*NUM_GNSS_CH - 1);
  RemoveColumns(S, NUM_GNSS_CH + observe_gnss_t, 2*NUM_GNSS_CH - 1);
  RemoveRows(S, observe_gnss_m, NUM_GNSS_CH - 1);
  RemoveColumns(S, observe_gnss_m, NUM_GNSS_CH - 1);
};

void PBD_dgps::ResizeMHt(Eigen::MatrixXd& MHt, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c)
{
  // ここも後ろから．
  RemoveRows(MHt, NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + observe_gnss_t, NUM_STATE_ALL - 1);
  RemoveRows(MHt, NUM_SINGLE_STATE + observe_gnss_m, NUM_SINGLE_STATE_ALL - 1);
  RemoveColumns(MHt, 2*NUM_GNSS_CH + observe_gnss_c, 3*NUM_GNSS_CH - 1);
  RemoveColumns(MHt, NUM_GNSS_CH + observe_gnss_t, 2*NUM_GNSS_CH - 1);
  RemoveColumns(MHt, observe_gnss_m, NUM_GNSS_CH - 1);
};

// FIXME: 再度この辺の観測行列にミスがないかを確認する．
void PBD_dgps::UpdateObservationsGRAPHIC(const int sat_id, EstimatedVariables& x_est, const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv, const int N_offset)
{
  // ここもLEO satが把握している誤差ありの情報．
  // std::find index of the observing gnss satellite
  PBD_GnssObservation& gnss_observation = gnss_observations_.at(sat_id);
  GnssObserveModel& observe_model = gnss_observed_models_.at(sat_id);

  const GnssObserveInfo& observe_info_ = gnss_observation.info_;
  const int index = conv_index_from_gnss_sat_id(observe_info_.now_observed_gnss_sat_id, gnss_sat_id);

  const GnssObservedValues& observed_val_ = gnss_observation.observed_values_;
  auto gnss_position = observed_val_.gnss_satellites_position.at(index);
  double gnss_clock = observed_val_.gnss_clock.at(index);
  // とりあえずL1を使う．
  double pseudo_range = observed_val_.L1_pseudo_range.at(index);
  const auto& carrier_phase = observed_val_.L1_carrier_phase.at(index);
  double carrier_phase_range = carrier_phase.first * x_est.lambda; // [m]で受け取る．

  libra::Vector<3> sat_position{ 0.0 };
  sat_position = ConvEigenVecToLibraVec(x_est.position);
  const double sat_clock = x_est.clock(0);
  // ここら辺はGnssObserveModel内に格納する．
  libra::Vector<3> receive_pos = ConvEigenVecToLibraVec(ConvCenterOfMassToReceivePos(x_est.position, antenna_pos_b_.at(sat_id), sat_info_.at(sat_id).dynamics));
  double geometric_range = gnss_observation.CalculateGeometricRange(receive_pos, gnss_position);
  observe_model.geometric_range.push_back(geometric_range);
  if (observe_model.geometric_range.size() != (index + 1)) abort();

  double pseudo_range_model = gnss_observation.CalculatePseudoRange(sat_position, gnss_position, sat_clock, gnss_clock);
  observe_model.pseudo_range_model.push_back(pseudo_range_model);
  if (observe_model.pseudo_range_model.size() != (index + 1)) abort();

  double carrier_phase_model = gnss_observation.CalculateCarrierPhase(sat_position, gnss_position, sat_clock, gnss_clock, x_est.ambiguity.N.at(index), L1_lambda);
  observe_model.carrier_phase_model.push_back(carrier_phase_model);
  if (observe_model.carrier_phase_model.size() != (index + 1)) abort();

  // GRAPHIC
  const int row_offset = N_offset + index;
  const int num_state_offset = sat_id*(NUM_SINGLE_STATE + N_offset);
  const int col_offset = num_state_offset + NUM_SINGLE_STATE + index;
  z(row_offset) = (pseudo_range + carrier_phase_range) / 2;
  h_x(row_offset) = (pseudo_range_model + carrier_phase_model) / 2;
  for (int j = 0; j < 3; ++j) {
    const int col_index = num_state_offset + j;
    // position
    H(row_offset, col_index) = (x_est.position(j) - gnss_position[j]) / geometric_range;
  }
  // clock
  H(row_offset, num_state_offset + 3) = 1.0;
  // Ambiguity
#ifndef N_DEBUG
  H(row_offset, col_offset) = 0.5 * x_est.lambda; // GRAPHIC N
#endif // N_DEBUG
  Rv(row_offset) = pow(pseudo_sigma / 2.0, 2.0) + pow(carrier_sigma / 2.0, 2.0);
};

void PBD_dgps::UpdateObservationsSDCP(const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
{
  // std::find for comon index
  const int common_index = conv_index_from_gnss_sat_id(common_observed_gnss_sat_id, gnss_sat_id);
  const int row_offset = visible_gnss_nums_.at(0) + visible_gnss_nums_.at(1) + common_index;
  PBD_GnssObservation& main_observation = gnss_observations_.at(0);
  PBD_GnssObservation& target_observation = gnss_observations_.at(1);

  const int main_index = conv_index_from_gnss_sat_id(main_observation.info_.now_observed_gnss_sat_id, gnss_sat_id);
  const int target_index = conv_index_from_gnss_sat_id(target_observation.info_.now_observed_gnss_sat_id, gnss_sat_id);
  const int col_offset_main = NUM_SINGLE_STATE + main_index;
  const int col_offset_target = NUM_STATE + visible_gnss_nums_.at(0) + target_index;
  const int num_main_state_all = NUM_SINGLE_STATE + visible_gnss_nums_.at(0);

  // とりあえずL1を使う．
  // main
  const auto& carrier_phase_main = main_observation.observed_values_.L1_carrier_phase.at(main_index);
  double carrier_phase_range_main = carrier_phase_main.first * x_est_main.lambda; // FIXME: x_est_をmaster情報にできるように修正する．
  // target
  const auto& carrier_phase_target = target_observation.observed_values_.L1_carrier_phase.at(target_index);
  double carrier_phase_range_target = carrier_phase_target.first * x_est_target.lambda;

  auto gnss_position = main_observation.observed_values_.gnss_satellites_position.at(main_index);

  GnssObserveModel& main_observe_model = gnss_observed_models_.at(0);
  GnssObserveModel& target_observe_model = gnss_observed_models_.at(1);

  // SDCP
  z(row_offset) = carrier_phase_range_target - carrier_phase_range_main;
  h_x(row_offset) = target_observe_model.carrier_phase_model.at(target_index) - main_observe_model.carrier_phase_model.at(main_index);
  for (int j = 0; j < 3; ++j) {
    // position
    // main
    H(row_offset, j) = (gnss_position[j] - x_est_main.position(j)) / main_observe_model.geometric_range.at(main_index);
    // target
    H(row_offset, num_main_state_all + j) = (x_est_target.position(j) - gnss_position[j]) / target_observe_model.geometric_range.at(target_index);
  }
  // clock
  H(row_offset, 3) = -1.0;
  H(row_offset, num_main_state_all + 3) = 1.0;
#ifndef N_DEBUG
  // Ambiguity
  H(row_offset, col_offset_main) = -1.0 * x_est_main.lambda; // SDCP N
  H(row_offset, col_offset_target) = 1.0 * x_est_target.lambda; // SDCP N
#endif // N_DEBUG

  Rv(row_offset) = pow(sqrt(2.0) * carrier_sigma, 2.0);
};

// この内部で観測量関連を　すべて更新する．TODO: 更新しなくていいものもまとめてる気がするから，したいものだけに限定した方がいいかも．
void PBD_dgps::UpdateObservations(Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
{
  // FIXME: このアクセス方法では汎用性は低い．
  const GnssObserveInfo& main_info_ = gnss_observations_.at(0).info_;
  const GnssObserveInfo& target_info_ = gnss_observations_.at(1).info_;
  std::vector<int> all_observed_gnss_ids = main_info_.now_observed_gnss_sat_id;
  all_observed_gnss_ids.insert(all_observed_gnss_ids.end(), target_info_.now_observed_gnss_sat_id.begin(), target_info_.now_observed_gnss_sat_id.end()); // concate
  sort(all_observed_gnss_ids.begin(), all_observed_gnss_ids.end());
  all_observed_gnss_ids.erase(unique(all_observed_gnss_ids.begin(), all_observed_gnss_ids.end()), all_observed_gnss_ids.end()); // unique

  for (const int& id : all_observed_gnss_ids)
  {
    auto it_main = std::find(main_info_.now_observed_gnss_sat_id.begin(), main_info_.now_observed_gnss_sat_id.end(), id);
    if (it_main != main_info_.now_observed_gnss_sat_id.end() && *it_main == id)
    {
      UpdateObservationsGRAPHIC(0, x_est_main, id, z, h_x, H, Rv, 0);
    }
    // if target
    auto it_target = std::find(target_info_.now_observed_gnss_sat_id.begin(), target_info_.now_observed_gnss_sat_id.end(), id);
    if (it_target != target_info_.now_observed_gnss_sat_id.end() && *it_target == id)
    {
      UpdateObservationsGRAPHIC(1, x_est_target, id, z, h_x, H, Rv, visible_gnss_nums_.at(0));
    }
    // if common
    auto it_common = std::find(common_observed_gnss_sat_id.begin(), common_observed_gnss_sat_id.end(), id);
    if (it_common != common_observed_gnss_sat_id.end() && *it_common == id)
    {
      UpdateObservationsSDCP(id, z, h_x, H, Rv);
    }
  }
}

static const int conv_index_from_gnss_sat_id(std::vector<int> observed_gnss_sat_id, const int gnss_sat_id)
{
  std::vector<int>::iterator itr = std::find(observed_gnss_sat_id.begin(), observed_gnss_sat_id.end(), gnss_sat_id);
  if (itr == observed_gnss_sat_id.end())
  {
    std::cout << "not found" << gnss_sat_id << std::endl;
    abort();
  }
  const int index = distance(observed_gnss_sat_id.begin(), itr);
  return index;
};

// sat_idはLEO衛星のこと
void PBD_dgps::FindCommonObservedGnss(const std::pair<int, int> sat_id_pair)
{
  // 初期化
  common_observed_gnss_sat_id.clear();
  common_observed_status.assign(num_of_gnss_satellites_, false);
  const int main_sat_id = sat_id_pair.first;
  const int target_sat_id = sat_id_pair.second;
  int common_index = 0;

  const PBD_GnssObservation& main_observation_ = gnss_observations_.at(main_sat_id);
  const PBD_GnssObservation& target_observation_ = gnss_observations_.at(target_sat_id);
  // ここはiterで取得でいいのかも？
  for (int i = 0; i < visible_gnss_nums_.at(0); ++i)
  {
    for (int j = 0; j < visible_gnss_nums_.at(1); ++j)
    {
      // ここでずれている？
      int gnss_sat_id = main_observation_.info_.now_observed_gnss_sat_id.at(i);
      // どっかで複数衛星にも拡張
      if (gnss_sat_id == target_observation_.info_.now_observed_gnss_sat_id.at(j))
      {
        common_observed_gnss_sat_id.push_back(gnss_sat_id);
        common_observed_status.at(gnss_sat_id) = true;
        // common_index_dict.at(gnss_sat_id) = common_index; // このindexはほんまに必要なのか？
        ++common_index;
        // pre_common_observing_ch = now_common_observing_ch; // ???
        //AllocateToCh(gnss_sat_id, now_common_observing_ch, common_free_ch);
        break;
      }
      // pre_common_observing_ch = now_common_observing_ch;
      // RemoveFromCh(gnss_sat_id, now_common_observing_ch, common_free_ch);
    }
  }
}

// 空のcommon freeにアクセスして死んだ．
void PBD_dgps::AllocateToCh(const int gnss_sat_id, std::map<const int, int>& observing_ch, std::vector<int>& free_ch)
{
  if (observing_ch.count(gnss_sat_id))
  {
    // 引き続き観測なので更新しない
  }
  else {
    int ch = free_ch.at(0);
    free_ch.erase(free_ch.begin());
    observing_ch[gnss_sat_id] = ch;
  }
};

void PBD_dgps::RemoveFromCh(const int gnss_sat_id, std::map<const int, int>& observing_ch, std::vector<int>& free_ch)
{
  if (observing_ch.count(gnss_sat_id))
  {
    int ch = observing_ch[gnss_sat_id];
    observing_ch.erase(gnss_sat_id);
    free_ch.push_back(ch);
  }
  else
  {
    // None
  }
};

Eigen::Vector3d PBD_dgps::ConvReceivePosToCenterOfMass(Eigen::Vector3d rec_pos, libra::Vector<3> antenna_pos_b, const Dynamics& dynamics)
{
  Eigen::Vector3d pos;
  Quaternion q_i2b = dynamics.GetQuaternion_i2b();
  libra::Vector<3> sat2ant_i = q_i2b.frame_conv_inv(antenna_pos_b);
  for (uint8_t i = 0; i < 3; i++) pos(i) = rec_pos(i) - sat2ant_i[i]; // 補正する．

  return pos;
}

Eigen::Vector3d PBD_dgps::ConvCenterOfMassToReceivePos(Eigen::Vector3d pos, libra::Vector<3> antenna_pos_b, const Dynamics& dynamics)
{
  Eigen::Vector3d receive_pos;
  Quaternion q_i2b = main_dynamics_.GetQuaternion_i2b();
  libra::Vector<3> sat2ant_i = q_i2b.frame_conv_inv(antenna_pos_b);
  for (uint8_t i = 0; i < 3; i++) receive_pos(i) = pos(i) + sat2ant_i[i]; // 補正する．

  return receive_pos;
}

void PBD_dgps::UpdateBiasForm(const int sat_id, EstimatedVariables& x_est, Eigen::MatrixXd& P, Eigen::MatrixXd& Q) // LEO衛星の数が増えたときは衛星ごとにこのクラスのインスタンスが生成される？ので一旦これで行く
{
  const PBD_GnssObservation& gnss_observation = gnss_observations_.at(sat_id);
  const GnssObserveInfo& observe_info_ = gnss_observation.info_;
  //観測する衛星同じだったら飛ばしていい
  if (CheckVectorEqual(gnss_observation.info_.pre_observed_gnss_sat_id, gnss_observation.info_.now_observed_gnss_sat_id))
  {
    return;
  }

  int n = visible_gnss_nums_.at(sat_id);
  int n_pre = pre_visible_gnss_nums_.at(sat_id);

  // index is the order in matrix
  int pre_index = 0;
  int now_index = 0;

  const std::vector<double> pre_estimated_bias = x_est.ambiguity.N;
  const Eigen::MatrixXd pre_P = P;
  const Eigen::MatrixXd pre_Q = Q;
  // reset
  const int num_state_all = NUM_SINGLE_STATE + n;
  const int pre_num_state_all = NUM_SINGLE_STATE + n_pre;
  P = Eigen::MatrixXd::Zero(num_state_all, num_state_all);
  Q = Eigen::MatrixXd::Zero(num_state_all, num_state_all);
  P.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE) = pre_P.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE);
  x_est.ambiguity.N.assign(NUM_GNSS_CH, 0); // resetする．
  Q.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE) = pre_Q.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE);

  // for debug
  std::vector<int> now_gnss_sat_ids = observe_info_.now_observed_gnss_sat_id;
  std::vector<int> pre_gnss_sat_ids = observe_info_.pre_observed_gnss_sat_id;

  for (int i = 0; i < now_gnss_sat_ids.size(); ++i)
  {
    const int gnss_sat_id = now_gnss_sat_ids.at(i);
    if (observe_info_.pre_observed_status.at(gnss_sat_id) == false && observe_info_.now_observed_status.at(gnss_sat_id) == false) continue; // もはやここは必要ない．
    // 見えなくなったとき
    else if (observe_info_.pre_observed_status.at(gnss_sat_id) == true && observe_info_.now_observed_status.at(gnss_sat_id) == false) // これもいらない．
    {
      // 何もせず飛ばす．
      if (pre_index != conv_index_from_gnss_sat_id(pre_gnss_sat_ids, gnss_sat_id)) abort();
      ++pre_index;
    }
    // else if (observe_info_.pre_observed_status.at(gnss_sat_id) == false && observe_info_.now_observed_status.at(gnss_sat_id) == true)
    else if (std::find(pre_gnss_sat_ids.begin(), pre_gnss_sat_ids.end(), gnss_sat_id) == pre_gnss_sat_ids.end() && observe_info_.now_observed_status.at(gnss_sat_id) == true)
    {
      if (now_index != conv_index_from_gnss_sat_id(now_gnss_sat_ids, gnss_sat_id)) abort();
      // std::normal_distribution<> N_dist(0.0, sigma_N_ini);

      Eigen::Vector3d x_est_rec = ConvCenterOfMassToReceivePos(x_est.position, antenna_pos_b_.at(sat_id), sat_info_.at(sat_id).dynamics);
      // 擬似距離観測量をそのまま使うバージョン
      double ionosphere_delay = gnss_observation.CalculateIonDelay(gnss_sat_id, ConvEigenVecToLibraVec(x_est_rec), L1_frequency); // 電離圏遅延量を既知とする．
      double observed_pseudo_range = gnss_observation.observed_values_.L1_pseudo_range.at(now_index) - ionosphere_delay;
      x_est.ambiguity.N.at(now_index) = (gnss_observation.observed_values_.L1_carrier_phase.at(now_index).first * x_est.lambda - observed_pseudo_range + ionosphere_delay) / x_est.lambda; // biasの初期値は搬送波位相距離と観測搬送波位相の差をとる．

      // モデルを使うバージョン
      // double pseudo_range_model = gnss_observation.CalculatePseudoRange(ConvEigenVecToLibraVec(x_est_rec), gnss_observation.observed_values_.gnss_satellites_position.at(now_index), x_est.clock(0), gnss_observation.observed_values_.gnss_clock.at(now_index)); // 本来はここに電離圏モデルを入れないとダメだが，フリーにしているので使いまわす．
      // x_est.ambiguity.N.at(now_index) = (gnss_observation.observed_values_.L1_carrier_phase.at(now_index).first * x_est.lambda - pseudo_range_model) / x_est.lambda; // biasの初期値は搬送波位相距離と観測搬送波位相の差をとる．

      sat_info_.at(sat_id).true_N(now_index) = - gnss_observation.l1_bias_.at(gnss_sat_id);
      int offset = NUM_SINGLE_STATE + now_index;

      // 行と列に関してもリセット
      P(offset, offset) = pow(sigma_N_ini, 2.0);
      Q(offset, offset) = pow(sigma_N_process, 2.0);
      ++now_index;
    }
    // 引き継ぐ
    // else if (observe_info_.pre_observed_status.at(gnss_sat_id) == true && observe_info_.now_observed_status.at(gnss_sat_id) == true)
    else if (std::find(pre_gnss_sat_ids.begin(), pre_gnss_sat_ids.end(), gnss_sat_id) != pre_gnss_sat_ids.end() && observe_info_.now_observed_status.at(gnss_sat_id) == true)
    {
      // if (pre_index != conv_index_from_gnss_sat_id(pre_gnss_sat_ids, gnss_sat_id)) abort();
      pre_index = conv_index_from_gnss_sat_id(pre_gnss_sat_ids, gnss_sat_id);
      if (now_index != conv_index_from_gnss_sat_id(now_gnss_sat_ids, gnss_sat_id)) abort();

      x_est.ambiguity.N.at(now_index) = pre_estimated_bias.at(pre_index);
      sat_info_.at(sat_id).true_N(now_index) = - gnss_observation.l1_bias_.at(gnss_sat_id); // 真値をとってくる．

      // 整数不定性以外との相関は引き継ぐ．
      P.block(0, NUM_SINGLE_STATE + now_index, NUM_SINGLE_STATE, 1) = pre_P.block(0, NUM_SINGLE_STATE + pre_index, NUM_SINGLE_STATE, 1);
      P.block(NUM_SINGLE_STATE + now_index, 0, 1, NUM_SINGLE_STATE) = pre_P.block(NUM_SINGLE_STATE + pre_index, 0, 1, NUM_SINGLE_STATE);

      // N間の関係は見ていた部分のみを引き継ぐ．
      for (int j = 0; j < now_gnss_sat_ids.size(); j++)
      {
        if (j == now_index) break; // 今の衛星以上の部分はまだ知らないのでここで終了．
        int other_gnss_sat_id = now_gnss_sat_ids.at(j);
        if (std::find(pre_gnss_sat_ids.begin(), pre_gnss_sat_ids.end(), other_gnss_sat_id) == pre_gnss_sat_ids.end()) continue; // 見えてなかった衛星に対してはスキップ
        int other_pre_index = conv_index_from_gnss_sat_id(pre_gnss_sat_ids, other_gnss_sat_id);

        P(NUM_SINGLE_STATE + j, NUM_SINGLE_STATE + now_index) = pre_P(NUM_SINGLE_STATE + other_pre_index, NUM_SINGLE_STATE + pre_index);
        P(NUM_SINGLE_STATE + now_index, NUM_SINGLE_STATE + j) = pre_P(NUM_SINGLE_STATE + pre_index, NUM_SINGLE_STATE + other_pre_index);
      }

      // 対角成分の引継ぎ
      P(NUM_SINGLE_STATE + now_index, NUM_SINGLE_STATE + now_index) = pre_P(NUM_SINGLE_STATE + pre_index, NUM_SINGLE_STATE + pre_index);
      Q(NUM_SINGLE_STATE + now_index, NUM_SINGLE_STATE + now_index) = pre_Q(NUM_SINGLE_STATE + pre_index, NUM_SINGLE_STATE + pre_index);
      ++now_index;
    }

    if (now_index >= NUM_GNSS_CH || pre_index >= NUM_GNSS_CH) break; // ch以上の受信は出来ない
  }

  // N_trueの更新が必要
  for (int i = visible_gnss_nums_.at(sat_id); i < NUM_GNSS_CH; i++)
  {
    sat_info_.at(sat_id).true_N(i) = 0;
  }
}

template <typename T> bool PBD_dgps::CheckVectorEqual(const std::vector<T>& a, const std::vector<T>& b)
{
    if(a.size() != b.size()) return false;
    for(int i = 0;i < a.size();++i){
        if(a.at(i) != b.at(i)) return false;
    }

    return true;
}

int PBD_dgps::SelectBaseGnssSatellite(Eigen::VectorXd N, Eigen::MatrixXd P_N)
{
  Eigen::VectorXd Var_N = P_N.diagonal();
  std::vector<double> var_N;
  for (int i = 0; i < Var_N.size(); ++i) var_N.push_back(Var_N(i));
  int gnss_sat_id;
  // 使うのはPの対角成分だけでいい．
  // TODO: 実装する．
  return gnss_sat_id;
}

void PBD_dgps::DynamicNoiseScaling(Eigen::MatrixXd Q_dash, Eigen::MatrixXd H)
{
  const int num_state_main = NUM_SINGLE_STATE + visible_gnss_nums_.at(0);
  const int num_state_all = num_state_main + NUM_SINGLE_STATE + visible_gnss_nums_.at(1);
  Eigen::MatrixXd Phi_all = Eigen::MatrixXd::Identity(num_state_all, num_state_all);

  Phi_all.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Phi_.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE);
  Phi_all.block(num_state_main, num_state_main, NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Phi_.bottomRightCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE);

  // clock
#ifdef CLOCK_IS_RANDOM_WALK
  Phi_all(3, 3) = 1.0; // random walkなのでΦは1
  Phi_all(num_state_main + 3, num_state_main + 3) = 1.0;
#endif // CLOCK_IS_RANDOM_WALK
#ifdef CLOCK_IS_WHITE_NOISE
  Phi_all(3, 3) = 0.0;
  Phi_all(num_state_main + 3, num_state_main + 3) = 0.0;
#endif // CLOCK_IS_WHITE_NOISE

  // a
#ifdef REDUCED_DYNAMIC
  Phi_all.block(7, 7, 3, 3) = CalculatePhi_a(observe_step_time);
  Phi_all.block(num_state_main + 7, num_state_main + 7, 3, 3) = CalculatePhi_a(observe_step_time);
#endif // REDUCED_DYNAMIC

  // 観測更新後のPに対して行う．
  Eigen::MatrixXd P_dash = Phi_all * P_ * Phi_all.transpose() + Q_dash;
  Eigen::MatrixXd P      = Phi_all * P_ * Phi_all.transpose() + Q_;
#define TRACE_SCALE_

#ifdef TRACE_SCALE_
  double beta_dash = (H * P_dash * H.transpose()).trace() / (H * P * H.transpose()).trace();
  Q_ = sqrt(beta_dash) * Q_; // Nは避けた方がよさそう．
#else
  // traceとらない方法.位置精度に依存する．
  Eigen::VectorXd diag_dash = P_dash.diagonal();
  Eigen::VectorXd diag = P.diagonal();
  for (uint8_t i = 0; i < diag.size(); i++)
  {
    if (diag(i) == 0) continue;
    Q_(i, i) = sqrt(diag_dash(i) / diag(i)) * Q_(i, i);
  }
#endif
}
