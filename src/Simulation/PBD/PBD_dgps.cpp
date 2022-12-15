#include <iomanip>
#include <set>
#include "PBD_dgps.h"
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
// #define TIME_UPDATE_DEBUG

#define PCO

#define LAMBDA_DEBUG

#undef cross

static const double PBD_DGPS_kConvNm2m = 1e-9;
static const int precision = 15; // position

// template使っていい感じにしたい．
// template <typename T> static void LogOutput(const )
static void LogOutput_(std::ofstream& ofs, const Eigen::MatrixXd& M, const int size, const int max_size);

PBD_dgps::PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const std::vector<PBD_Sat*> spacecrafts, PBD_GeoPotential* geop) :mt(42), step_time(sim_time_.GetStepSec()), ofs("result_new.csv"), num_of_gnss_satellites_(gnss_satellites_.GetNumOfSatellites()), geo_potential_(geop)
{
  //初期化
  x_est_main.position = Eigen::VectorXd::Zero(3);
  const Dynamics& main_dynamics = spacecrafts.at(0)->GetDynamics();
  libra::Vector<3> position_main = main_dynamics.GetPosition_i();
  for (int i = 0; i < 3; ++i) x_est_main.position(i) = position_main[i];
  x_est_main.clock = Eigen::VectorXd::Zero(1);
  libra::Vector<3> velocity_main = main_dynamics.GetOrbit().GetSatVelocity_i();
  for (int i = 0; i < 3; ++i) x_est_main.velocity(i) = velocity_main[i];
  x_est_main.acceleration = Eigen::VectorXd::Zero(3);
  x_est_main.acc_dist = Eigen::VectorXd::Zero(3);
  InitAmbiguity(x_est_main);
  const libra::Vector<3> pco_true = spacecrafts.at(0)->pbd_components_->GetGNSSReceiver()->GetPCO_mm();
  x_est_main.pcc = new PhaseCenterCorrection(pco_true, 5, 5);

  x_est_target.position = Eigen::VectorXd::Zero(3);
  const Dynamics& target_dynamics = spacecrafts.at(1)->GetDynamics();
  libra::Vector<3> position_target = target_dynamics.GetPosition_i();
  for (int i = 0; i < 3; ++i) x_est_target.position(i) = position_target[i];
  x_est_target.clock = Eigen::VectorXd::Zero(1);
  libra::Vector<3> velocity_target = target_dynamics.GetOrbit().GetSatVelocity_i();
  for (int i = 0; i < 3; ++i) x_est_target.velocity(i) = velocity_target[i];
  x_est_target.acceleration = Eigen::VectorXd::Zero(3);
  x_est_target.acc_dist = Eigen::VectorXd::Zero(3);
  InitAmbiguity(x_est_target);
  libra::Vector<3> initial_pco(0); initial_pco[2] = 120.0;
  x_est_target.pcc = new PhaseCenterCorrection(initial_pco, 5, 5);
  // x_est_target.pcc = new PhaseCenterCorrection(pco_true, 5, 5);

  // x_est_.push_back(x_est_main);
  // x_est_.push_back(x_est_target);
  sat_info_.push_back({main_dynamics, x_est_main, Eigen::VectorXd::Zero(NUM_GNSS_CH), spacecrafts.at(0)->pbd_components_});
  sat_info_.push_back({target_dynamics, x_est_target, Eigen::VectorXd::Zero(NUM_GNSS_CH), spacecrafts.at(1)->pbd_components_});

  gnss_observations_.push_back(*(spacecrafts.at(0)->gnss_observation_));
  gnss_observations_.push_back(*(spacecrafts.at(1)->gnss_observation_));

  GnssObserveModel main_model{};
  GnssObserveModel target_model{};
  gnss_observed_models_.push_back(main_model);
  gnss_observed_models_.push_back(target_model);

  antenna_pos_b_.push_back(gnss_observations_.at(0).GetAntennaPosition());
  antenna_pos_b_.push_back(gnss_observations_.at(1).GetAntennaPosition());

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
  libra::Vector<3> alignment_err_main = gnss_observations_.at(0).GetAntennaAlignmentError();
  libra::Vector<3> alignment_err_target = gnss_observations_.at(1).GetAntennaAlignmentError();
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
    const Dynamics& main_dynamics = sat_info_.at(0).dynamics;
    libra::Vector<3> sat_position_main = main_dynamics.GetPosition_i();
    libra::Vector<3> sat_velocity_main = main_dynamics.GetOrbit().GetSatVelocity_i();
    const Dynamics& target_dynamics = sat_info_.at(1).dynamics;
    libra::Vector<3> sat_position_target = target_dynamics.GetPosition_i();
    libra::Vector<3> sat_velocity_target = target_dynamics.GetOrbit().GetSatVelocity_i();

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

    // PCO
    const libra::Vector<3> pco_main = x_est_main.pcc->GetPCO_mm();
    for (int i = 0; i < 3; i++) ofs << std::fixed << std::setprecision(precision) << pco_main[i] << ",";
    const libra::Vector<3> pco_target = x_est_target.pcc->GetPCO_mm();
    for (int i = 0; i < 3; i++) ofs << std::fixed << std::setprecision(precision) << pco_target[i] << ",";

    // is_fixedフラグ
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i < visible_gnss_nums_.at(0)) ofs << std::fixed << std::setprecision(precision) << x_est_main.ambiguity.is_fixed.at(i) << ",";
      else ofs << false << ",";
    }
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i < visible_gnss_nums_.at(1)) ofs << std::fixed << std::setprecision(precision) << x_est_target.ambiguity.is_fixed.at(i) << ","; // N_est
      else ofs << false << ",";
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
  // clock部分
  Phi_(3, 3) = 0.0;
  Phi_(NUM_SINGLE_STATE + 3, NUM_SINGLE_STATE + 3) = 0.0;
  // scaling matrixの部分
  Phi_.block(7, 7, 3, 3) = Eigen::Matrix3d::Zero();
  Phi_.block(NUM_SINGLE_STATE + 7, NUM_SINGLE_STATE + 7, 3, 3) = Eigen::Matrix3d::Zero();
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
  Eigen::MatrixXd Phi_main = Phi_.topLeftCorner(7, NUM_SINGLE_STATE);
  Eigen::MatrixXd Phi_target = Phi_.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, 7, NUM_SINGLE_STATE);

  RK4(x_est_main.position, x_est_main.velocity, x_est_main.acceleration, x_est_main.acc_dist, Phi_main);
  Phi_.topLeftCorner(7, NUM_SINGLE_STATE) = Phi_main;

  RK4(x_est_target.position, x_est_target.velocity, x_est_target.acceleration, x_est_target.acc_dist, Phi_target);
  Phi_.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, 7, NUM_SINGLE_STATE) = Phi_target;

#ifdef REDUCED_DYNAMIC
  // acceleration ここでは q = sigma_acc_process
  Eigen::Matrix3d Phi_a = CalculatePhi_a(observe_step_time);
  double phi = Phi_a(0, 0);
  // double phi_t = Phi_a(1, 1);
  // double phi_n = Phi_a(2, 2);
  // x_est_main.acceleration = Phi_a * x_est_main.acceleration;
  // x_est_target.acceleration = Phi_a * x_est_target.acceleration;
  Phi_.block(7, 7, 3, 3) = Phi_a;
  Phi_.block(NUM_SINGLE_STATE + 7, NUM_SINGLE_STATE + 7, 3, 3) = Phi_a;

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
  AddGeoPotentialDisturbance(position, acc_dist);
  Eigen::Vector3d acc_dist_cpy = acc_dist;
  Eigen::Vector3d k0 = PositionDifferential(velocity);
  Eigen::Vector3d l0 = VelocityDifferential(position, velocity, acceleration, acc_dist);
  Eigen::Vector3d m0 = -acceleration / tau_a;
  Eigen::Matrix<double, 7, NUM_SINGLE_STATE> n0 = CalculateJacobian(position, velocity)*Phi;
  // (a, p)
#ifdef REDUCED_DYNAMIC
  Eigen::Matrix3d rtn2eci = TransRTN2ECI(position, velocity);
  n0.block(4, 7, 3, 3) += PBD_DGPS_kConvNm2m*rtn2eci;
#endif // REDUCED_DYNAMIC

  Eigen::Vector3d tmp_position = position + k0 * step_time * 0.5;
  Eigen::Vector3d tmp_velocity = velocity + l0 * step_time * 0.5;
  Eigen::Vector3d tmp_acceleration = acceleration + m0 * step_time * 0.5;
  Eigen::Matrix<double, 7, NUM_SINGLE_STATE> tmp_Phi = Phi + n0 * step_time * 0.5;
  acc_dist = acc_dist_cpy;
  Eigen::Vector3d k1 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l1 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Vector3d m1 = -tmp_acceleration / tau_a;
  Eigen::Matrix<double, 7, NUM_SINGLE_STATE> n1 = CalculateJacobian(tmp_position, tmp_velocity)*tmp_Phi;
#ifdef REDUCED_DYNAMIC
  n1.block(4, 7, 3, 3) += PBD_DGPS_kConvNm2m*rtn2eci;
#endif // REDUCED_DYNAMIC

  tmp_position = position + k1 * step_time * 0.5;
  tmp_velocity = velocity + l1 * step_time * 0.5;
  tmp_acceleration = acceleration + m1 * step_time * 0.5;
  tmp_Phi = Phi + n1 * step_time * 0.5;
  acc_dist = acc_dist_cpy;
  Eigen::Vector3d k2 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l2 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Vector3d m2 = -tmp_acceleration / tau_a;
  Eigen::Matrix<double, 7, NUM_SINGLE_STATE> n2 = CalculateJacobian(tmp_position, tmp_velocity)*tmp_Phi;
#ifdef REDUCED_DYNAMIC
  n2.block(4, 7, 3, 3) += PBD_DGPS_kConvNm2m*rtn2eci;
#endif // REDUCED_DYNAMIC

  tmp_position = position + k2 * step_time;
  tmp_velocity = velocity + l2 * step_time;
  tmp_acceleration = acceleration + m2 * step_time;
  tmp_Phi = Phi + n2 * step_time;
  acc_dist = acc_dist_cpy;
  Eigen::Vector3d k3 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l3 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Vector3d m3 = -tmp_acceleration / tau_a;
  Eigen::Matrix<double, 7, NUM_SINGLE_STATE> n3 = CalculateJacobian(tmp_position, tmp_velocity)*tmp_Phi;
#ifdef REDUCED_DYNAMIC
  n3.block(4, 7, 3, 3) += PBD_DGPS_kConvNm2m*rtn2eci;
#endif // REDUCED_DYNAMIC

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
  // AddGeoPotentialDisturbance(position, acc_dist); // coreに合わせるために一旦コメントアウトする．
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


  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(7, 7);
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
  Eigen::VectorXd pre_Rv = Eigen::VectorXd::Zero(num_observables); // 観測誤差共分散．
  for (int i = 0; i < num_observables; i++) pre_Rv(i) = R_(i, i); // このままだと観測衛星数が変わる場合に対応できない．
  UpdateObservations(z, h_x, H, pre_Rv);

  // ここにfixしたNの部分が含まれているから不安定になる？
  Eigen::MatrixXd hmh = H * P_ * H.transpose();

  Eigen::MatrixXd tmp = R_ + hmh; // (observation_num, observation_num)
  //カルマンゲイン
  Eigen::MatrixXd K = CalculateK(H, tmp);

  Eigen::VectorXd x_predict(num_state_all);
  x_predict <<  x_est_main.position, x_est_main.clock, x_est_main.velocity,
#ifdef REDUCED_DYNAMIC
                x_est_main.acceleration,
#endif
                ConvStdVecToEigenVec(x_est_main.ambiguity.N),
                x_est_target.position, x_est_target.clock, x_est_target.velocity,
#ifdef REDUCED_DYNAMIC
                x_est_target.acceleration,
#endif
                ConvStdVecToEigenVec(x_est_target.ambiguity.N);

  // innovationの記号を何にするかは要検討
  Eigen::VectorXd E_pre = z - h_x;

  Eigen::VectorXd x_ant_predict = x_predict;
#ifndef RANGE_OBSERVE_DEBUG
  // まずアンテナ位置に変換
  x_ant_predict.topRows(3) = ConvCenterOfMassToReceivePos(x_predict.topRows(3), antenna_pos_b_.at(0), sat_info_.at(0).dynamics);
  x_ant_predict.block(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), 0, 3, 1) = ConvCenterOfMassToReceivePos(x_predict.block(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), 0, 3, 1), antenna_pos_b_.at(1), sat_info_.at(1).dynamics);
#endif // RANGE_OBSERVE_DEBUG
  Eigen::VectorXd x_update = x_ant_predict + K * E_pre;
  // Eigen::VectorXd x_update = x_predict + K * E_pre;
#ifndef RANGE_OBSERVE_DEBUG
  // 重心位置に戻す．
  x_update.topRows(3) = ConvReceivePosToCenterOfMass(x_update.topRows(3), antenna_pos_b_.at(0), sat_info_.at(0).dynamics);
  x_update.block(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), 0, 3, 1) = ConvReceivePosToCenterOfMass(x_update.block(NUM_SINGLE_STATE + visible_gnss_nums_.at(0), 0, 3, 1), antenna_pos_b_.at(1), sat_info_.at(1).dynamics);
#endif // RANGE_OBSERVE_DEBUG

  //更新
  GetStateFromVector(num_main_state_all, x_update); // N以外を代入

  for (int i = 0; i < NUM_GNSS_CH; ++i)
  {
    if (i < n_main && !x_est_main.ambiguity.is_fixed.at(i)) x_est_main.ambiguity.N.at(i) = x_update(NUM_SINGLE_STATE + i);
    if (i < n_target && !x_est_target.ambiguity.is_fixed.at(i)) x_est_target.ambiguity.N.at(i) = x_update(num_main_state_all + NUM_SINGLE_STATE + i);
  }

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(num_state_all, num_state_all);
  tmp = (I - K*H);
  P_ = tmp*P_*tmp.transpose() + K*R_*K.transpose(); // ここでNの部分にも入ってくる？
  // P_ = tmp*P_;

#ifdef AKF
  // residual
  // GRAPHIC*2 + SDCPにする．
  Eigen::VectorXd h_x_post = Eigen::VectorXd::Zero(num_observables); // 観測モデル行列

  for (int i = 0; i < num_observables; i++) pre_Rv(i) = R_(i, i); // このままだと観測衛星数が変わる場合に対応できない．
  UpdateObservations(z, h_x_post, H, pre_Rv);

  Eigen::VectorXd E_post = z - h_x_post;
  // R_ = alpha * R_ + (1 - alpha) * (E_post*E_post.transpose() + H*P_*H.transpose()); // residual based R_-adaptation
  R_ = alpha * R_ + (1 - alpha) * E_post*E_post.transpose(); // residual based R_-adaptation
  Eigen::MatrixXd Q_dash = alpha * Q_ + (1 - alpha) * K * E_pre * (K * E_pre).transpose(); // Innovation based Q_-adaptation

  DynamicNoiseScaling(Q_dash, H);
  // これしてるから，絶対軌道精度の影響（つまり残差の大きさ）ではなくて，収束してしまう？
#endif // AKF

  // IAR
#ifdef LAMBDA_DEBUG
  bool lambda_result = false;
  Eigen::MatrixXd P_N_main = P_.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, n_main, n_main);
  Eigen::MatrixXd P_N_target = P_.block(NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE, NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE, n_target, n_target);
  if (P_N_main.maxCoeff() < 1.0 && P_N_target.maxCoeff() < 1.0)
  {
    std::vector<Eigen::MatrixXd> vec_P_N;
    vec_P_N.push_back(P_N_main);
    vec_P_N.push_back(P_N_target);
    std::vector<Ambiguity> vec_N_cpy{x_est_main.ambiguity, x_est_target.ambiguity}; // 失敗した時に戻す用
    std::vector<Ambiguity*> vec_N{&(x_est_main.ambiguity), &(x_est_target.ambiguity)};
    std::vector<std::vector<int>> observed_gnss_ids{gnss_observations_.at(0).info_.now_observed_gnss_sat_id, gnss_observations_.at(1).info_.now_observed_gnss_sat_id, common_observed_gnss_sat_id}; // ちょっと煩雑やけどまあ．．

    PBD_Lambda lambda(vec_N, vec_P_N, observed_gnss_ids, NUM_SINGLE_STATE);
    lambda_result = lambda.Solve(METHOD::PAR_ILS);
    if (lambda_result)
    {
      // 解けたNをベースに他状態量の更新を実施する．
      const Eigen::MatrixXd& M2M1 = lambda.M2M1_;
      Eigen::MatrixXd P_fixed = M2M1 * P_ * M2M1.transpose();
      Eigen::VectorXd b_hat(NUM_STATE);
      b_hat << x_update.topRows(NUM_SINGLE_STATE), x_update.block(NUM_SINGLE_STATE + n_main, 0, NUM_SINGLE_STATE, 1);
      Eigen::VectorXd b_fixed = b_hat - P_fixed.block(0, NUM_STATE, NUM_STATE, n_common - 1) * P_fixed.bottomRightCorner(n_common - 1, n_common - 1).inverse() * (lambda.N_hat_ - lambda.N_);
      Eigen::VectorXd x_fixed_all(NUM_STATE + n_main + n_target);
      x_fixed_all << b_fixed.topRows(NUM_SINGLE_STATE), ConvStdVecToEigenVec(vec_N.at(0)->N), b_fixed.bottomRows(NUM_SINGLE_STATE), ConvStdVecToEigenVec(vec_N.at(1)->N);
      // 代入
      GetStateFromVector(NUM_SINGLE_STATE + n_main, x_fixed_all); // Nはポインタで渡しているので必要なし．

      // 残差確認
      Eigen::VectorXd h_x_post_lambda = Eigen::VectorXd::Zero(num_observables); // 観測モデル行列
      UpdateObservations(z, h_x_post_lambda, H, pre_Rv); // 更新するのはSDCPだけでいいので無駄．
      Eigen::MatrixXd E_post_lambda = z - h_x_post_lambda;
      if ((abs(E_post_lambda.bottomRows(n_common).array()) > L1_lambda).any()) // SDCPだけを確認する．
      {
        std::cout << "false integer!" << std::endl;
        GetStateFromVector(NUM_SINGLE_STATE + n_main, x_update); // 戻す．
        x_est_main.ambiguity = vec_N_cpy.at(0); // falseにセットしないとダメそう．
        x_est_target.ambiguity = vec_N_cpy.at(1);
        lambda_result = false;
      }
      else
      {
        // 解けたらNの分散を落とす．
        for (int i = 0; i < n_main; i++)
        {
          if (!x_est_main.ambiguity.is_fixed.at(i)) continue;
          // P_(NUM_SINGLE_STATE + i, NUM_SINGLE_STATE + i) = 0;
          P_.col(NUM_SINGLE_STATE + i) = Eigen::MatrixXd::Zero(num_state_all, 1);
          P_.row(NUM_SINGLE_STATE + i) = Eigen::MatrixXd::Zero(1, num_state_all);
        }
        for (int i = 0; i < n_target; i++)
        {
          if (!x_est_target.ambiguity.is_fixed.at(i)) continue;
          // P_(num_main_state_all + NUM_SINGLE_STATE + i, num_main_state_all + NUM_SINGLE_STATE + i) = 0;
          P_.col(num_main_state_all + NUM_SINGLE_STATE + i) = Eigen::MatrixXd::Zero(num_state_all, 1);
          P_.row(num_main_state_all + NUM_SINGLE_STATE + i) = Eigen::MatrixXd::Zero(1, num_state_all);
        }
        // process noiseも同様．
        // Q_.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, n_main, n_main) = Eigen::MatrixXd::Zero(n_main, n_main);
        // Q_.block(num_main_state_all + NUM_SINGLE_STATE, num_main_state_all + NUM_SINGLE_STATE, n_target, n_target) = Eigen::MatrixXd::Zero(n_target, n_target);
      }

// PCOの推定．
#ifdef PCO
      // 全部fixしているときに限定する．
      if (std::count(x_est_main.ambiguity.is_fixed.begin(), x_est_main.ambiguity.is_fixed.end(), true) == n_main &&
          std::count(x_est_target.ambiguity.is_fixed.begin(), x_est_target.ambiguity.is_fixed.end(), true) == n_target)
      {
        // 推定完了したら実施しない．
        if (!x_est_target.pcc->GetPcoFixed()) EstimateDeltaPCO(ConvEigenVecToStdVec(z.bottomRows(n_common)));
      }
#endif // PCO

    }
  }
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
  Eigen::MatrixXd PHt = P_ * H.transpose();

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

// ここら辺template使って整理したい．
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

void PBD_dgps::UpdateObservationsGRAPHIC(const int sat_id, EstimatedVariables& x_est, const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv, const int N_offset)
{
  // ここもLEO satが把握している誤差ありの情報．
  // std::find index of the observing gnss satellite
  PBD_GnssObservation& gnss_observation = gnss_observations_.at(sat_id);
  GnssObserveModel& observe_model = gnss_observed_models_.at(sat_id);

  const GnssObserveInfo& observe_info_ = gnss_observation.info_;
  const int index = GetIndexOfStdVector<int>(observe_info_.now_observed_gnss_sat_id, gnss_sat_id);

  const GnssObservedValues& observed_val_ = gnss_observation.observed_values_;
  // とりあえずL1を使う．
  double pseudo_range = observed_val_.L1_pseudo_range.at(index);
  const auto& carrier_phase = observed_val_.L1_carrier_phase.at(index);
  double carrier_phase_range = carrier_phase.first * x_est.lambda; // [m]で受け取る．

  libra::Vector<3> sat_position = ConvEigenVecToLibraVec<3>(x_est.position);
  const double sat_clock = x_est.clock(0);
  // ここら辺はGnssObserveModel内に格納する．
  // libra::Vector<3> receive_pos = ConvEigenVecToLibraVec<3>(ConvCenterOfMassToReceivePos(x_est.position, antenna_pos_b_.at(sat_id), sat_info_.at(sat_id).dynamics));
  double geometric_range = gnss_observation.CalculateGeometricRange(gnss_sat_id, sat_position);
  observe_model.geometric_range.push_back(geometric_range);
  if (observe_model.geometric_range.size() != (index + 1)) abort();

  double pseudo_range_model = gnss_observation.CalculatePseudoRange(gnss_sat_id, sat_position, sat_clock);
  observe_model.pseudo_range_model.push_back(pseudo_range_model);
  if (observe_model.pseudo_range_model.size() != (index + 1)) abort();

  const double pcc = x_est.pcc->GetPCC_m(gnss_observation.GetGnssAzimuthDeg(index), gnss_observation.GetGnssElevationDeg(index));
  double carrier_phase_model = gnss_observation.CalculateCarrierPhase(gnss_sat_id, sat_position, sat_clock, x_est.ambiguity.N.at(index), L1_lambda, pcc);
  observe_model.carrier_phase_model.push_back(carrier_phase_model);
  if (observe_model.carrier_phase_model.size() != (index + 1)) abort();

  // GRAPHIC
  const int row_offset = N_offset + index;
  const int num_state_offset = sat_id*(NUM_SINGLE_STATE + N_offset);
  const int col_offset = num_state_offset + NUM_SINGLE_STATE + index;
  z(row_offset) = (pseudo_range + carrier_phase_range) / 2;
  h_x(row_offset) = (pseudo_range_model + carrier_phase_model) / 2;
  const libra::Vector<3> los_vec = gnss_observation.GetGnssDirection_i(index);
  for (int j = 0; j < 3; ++j) {
    const int col_index = num_state_offset + j;
    // position
    H(row_offset, col_index) = - los_vec[j];
  }
  // clock
  H(row_offset, num_state_offset + 3) = 1.0;
  // Ambiguity
#ifndef N_DEBUG
  H(row_offset, col_offset) = 0.5 * x_est.lambda; // GRAPHIC N
#endif // N_DEBUG
};

void PBD_dgps::UpdateObservationsSDCP(const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
{
  const int common_index = GetIndexOfStdVector<int>(common_observed_gnss_sat_id, gnss_sat_id);
  const int row_offset = visible_gnss_nums_.at(0) + visible_gnss_nums_.at(1) + common_index;
  PBD_GnssObservation& main_observation = gnss_observations_.at(0);
  PBD_GnssObservation& target_observation = gnss_observations_.at(1);

  const int main_index = GetIndexOfStdVector<int>(main_observation.info_.now_observed_gnss_sat_id, gnss_sat_id);
  const int target_index = GetIndexOfStdVector<int>(target_observation.info_.now_observed_gnss_sat_id, gnss_sat_id);
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

  // SDCP
  z(row_offset) = carrier_phase_range_target - carrier_phase_range_main;
  h_x(row_offset) = gnss_observed_models_.at(1).carrier_phase_model.at(target_index) - gnss_observed_models_.at(0).carrier_phase_model.at(main_index);
  const libra::Vector<3> los_vec_main = main_observation.GetGnssDirection_i(main_index);
  const libra::Vector<3> los_vec_target = target_observation.GetGnssDirection_i(target_index);
  for (int j = 0; j < 3; ++j) {
    // position
    // main
    H(row_offset, j) = los_vec_main[j];
    // target
    H(row_offset, num_main_state_all + j) = - los_vec_target[j];
  }
  // clock
  H(row_offset, 3) = -1.0;
  H(row_offset, num_main_state_all + 3) = 1.0;
#ifndef N_DEBUG
  // Ambiguity
  H(row_offset, col_offset_main) = -1.0 * x_est_main.lambda; // SDCP N
  H(row_offset, col_offset_target) = 1.0 * x_est_target.lambda; // SDCP N
#endif // N_DEBUG
};

// この内部で観測量関連を　すべて更新する．TODO: 更新しなくていいものもまとめてる気がするから，したいものだけに限定した方がいいかも．
void PBD_dgps::UpdateObservations(Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
{
  // 追加する前に前のものをクリア．
  ClearGnssObserveModels(gnss_observed_models_.at(0));
  ClearGnssObserveModels(gnss_observed_models_.at(1));
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
      AdjustReceiveCovariance(main_info_.now_observed_gnss_sat_id, main_info_.pre_observed_gnss_sat_id, id, 0, Rv);
    }
    // if target
    auto it_target = std::find(target_info_.now_observed_gnss_sat_id.begin(), target_info_.now_observed_gnss_sat_id.end(), id);
    if (it_target != target_info_.now_observed_gnss_sat_id.end() && *it_target == id)
    {
      UpdateObservationsGRAPHIC(1, x_est_target, id, z, h_x, H, Rv, visible_gnss_nums_.at(0));
      AdjustReceiveCovariance(target_info_.now_observed_gnss_sat_id, target_info_.pre_observed_gnss_sat_id, id, visible_gnss_nums_.at(0), Rv);
    }
    // if common
    auto it_common = std::find(common_observed_gnss_sat_id.begin(), common_observed_gnss_sat_id.end(), id);
    if (it_common != common_observed_gnss_sat_id.end() && *it_common == id)
    {
      UpdateObservationsSDCP(id, z, h_x, H, Rv);
      AdjustReceiveCovariance(common_observed_gnss_sat_id, main_info_.pre_observed_gnss_sat_id, id, visible_gnss_nums_.at(0) + visible_gnss_nums_.at(1), Rv);
    }
  }
}

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
        break;
      }
      // pre_common_observing_ch = now_common_observing_ch;
    }
  }
}

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
  Quaternion q_i2b = dynamics.GetQuaternion_i2b();
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

  const Ambiguity pre_estimated_ambiguity = x_est.ambiguity;
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
      if (pre_index != GetIndexOfStdVector<int>(pre_gnss_sat_ids, gnss_sat_id)) abort();
      ++pre_index;
    }
    // else if (observe_info_.pre_observed_status.at(gnss_sat_id) == false && observe_info_.now_observed_status.at(gnss_sat_id) == true)
    else if (std::find(pre_gnss_sat_ids.begin(), pre_gnss_sat_ids.end(), gnss_sat_id) == pre_gnss_sat_ids.end() && observe_info_.now_observed_status.at(gnss_sat_id) == true)
    {
      if (now_index != GetIndexOfStdVector<int>(now_gnss_sat_ids, gnss_sat_id)) abort();
      // std::normal_distribution<> N_dist(0.0, sigma_N_ini);

      Eigen::Vector3d x_est_rec = ConvCenterOfMassToReceivePos(x_est.position, antenna_pos_b_.at(sat_id), sat_info_.at(sat_id).dynamics);
      double ionosphere_delay = gnss_observation.CalculateIonDelay(gnss_sat_id, ConvEigenVecToLibraVec<3>(x_est_rec), L1_frequency); // 電離圏遅延量を既知とする．

      // 擬似距離観測量をそのまま使うバージョン
        double observed_pseudo_range = gnss_observation.observed_values_.L1_pseudo_range.at(now_index) - ionosphere_delay;
        x_est.ambiguity.N.at(now_index) = (gnss_observation.observed_values_.L1_carrier_phase.at(now_index).first * x_est.lambda - observed_pseudo_range + ionosphere_delay) / x_est.lambda; // biasの初期値は搬送波位相距離と観測搬送波位相の差をとる．
        x_est.ambiguity.gnss_sat_id.at(now_index) = gnss_sat_id;
        x_est.ambiguity.is_fixed.at(now_index) = false;

      // モデルを使うバージョン (こっちだと絶対精度が悪くなってしまう．初期誤差の影響を受けすぎる感じ．)
      // double pseudo_range_model = gnss_observation.CalculatePseudoRange(ConvEigenVecToLibraVec<3>(x_est_rec), gnss_observation.observed_values_.gnss_satellites_position.at(now_index), x_est.clock(0), gnss_observation.observed_values_.gnss_clock.at(now_index)); // 本来はここに電離圏モデルを入れないとダメだが，フリーにしているので使いまわす．
      // x_est.ambiguity.N.at(now_index) = (gnss_observation.observed_values_.L1_carrier_phase.at(now_index).first * x_est.lambda - pseudo_range_model) / x_est.lambda; // biasの初期値は搬送波位相距離と観測搬送波位相の差をとる．

      sat_info_.at(sat_id).true_N(now_index) = - gnss_observation.l1_bias_.at(gnss_sat_id);
      int offset = NUM_SINGLE_STATE + now_index;

      // 対角成分だけ初期化．
      P(offset, offset) = pow(sigma_N_ini, 2.0);
      Q(offset, offset) = pow(sigma_N_process, 2.0);
      ++now_index;
    }
    // 引き継ぐ
    // else if (observe_info_.pre_observed_status.at(gnss_sat_id) == true && observe_info_.now_observed_status.at(gnss_sat_id) == true)
    else if (std::find(pre_gnss_sat_ids.begin(), pre_gnss_sat_ids.end(), gnss_sat_id) != pre_gnss_sat_ids.end() && observe_info_.now_observed_status.at(gnss_sat_id) == true)
    {
      // if (pre_index != GetIndexOfStdVector<int>(pre_gnss_sat_ids, gnss_sat_id)) abort();
      pre_index = GetIndexOfStdVector<int>(pre_gnss_sat_ids, gnss_sat_id);
      if (now_index != GetIndexOfStdVector<int>(now_gnss_sat_ids, gnss_sat_id)) abort();

      // ここを見ると一気に代入vができないのでデータ構造変えた方がよさそう．
      x_est.ambiguity.N.at(now_index) = pre_estimated_ambiguity.N.at(pre_index);
      x_est.ambiguity.gnss_sat_id.at(now_index) = gnss_sat_id;
      x_est.ambiguity.is_fixed.at(now_index) = pre_estimated_ambiguity.is_fixed.at(pre_index);

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
        int other_pre_index = GetIndexOfStdVector<int>(pre_gnss_sat_ids, other_gnss_sat_id);

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

// もう少しうまくやりたい．
void PBD_dgps::AdjustReceiveCovariance(const std::vector<int>& now_gnss_sat_ids, const std::vector<int>& pre_gnss_sat_ids, const int gnss_sat_id, const int base_offset, const Eigen::VectorXd& pre_Rv)
{
  const int now_index = GetIndexOfStdVector(now_gnss_sat_ids, gnss_sat_id);
  const int offset = base_offset + now_index;
  // 対角成分以外は0とする．
  R_.row(offset) = Eigen::MatrixXd::Zero(1, R_.cols());
  R_.col(offset) = Eigen::MatrixXd::Zero(R_.rows(), 1);
  if (std::find(pre_gnss_sat_ids.begin(), pre_gnss_sat_ids.end(), gnss_sat_id) != pre_gnss_sat_ids.end())
  {
    const int pre_index = GetIndexOfStdVector(pre_gnss_sat_ids, gnss_sat_id);
    R_(offset, offset) = pre_Rv(base_offset + pre_index); // 対角以外はリセットした方がいい？
  }
  else
  {
    // GRAPHIC
    if (base_offset <= NUM_GNSS_CH) R_(offset, offset) = pow(pseudo_sigma, 2.0) + pow(carrier_sigma, 2.0);
    // SDCP
    else R_(offset, offset) = 2.0 * pow(carrier_sigma, 2.0);
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

// クラスが肥大化してしまっているので分割したい．
void PBD_dgps::EstimateDeltaPCO(const std::vector<double> sdcp_vec)
{
  const PBD_GnssObservation& main_observation = gnss_observations_.at(0);
  const std::vector<GnssInfo> main_vec_gnssinfo = main_observation.GetReceiver()->GetGnssInfoVec();

  const RelativePositionSensor* rel_sensor = sat_info_.at(0).components->GetRelativePositionSensor();

  // 精密な相対位置を取得
  libra::Vector<3> relative_position_rtn = rel_sensor->GetMeasuredTargetPosition_rtn_m();
  Eigen::Vector3d relative_position_eci = TransRTN2ECI(x_est_main.position, x_est_main.velocity) * ConvLibraVecToEigenVec<3>(relative_position_rtn);

  x_est_target.position = x_est_main.position + relative_position_eci; // 正確な相対位置に更新．<-絶対位置がずれていても以下で差分をとるので良い．

  double max_elevation = 0;
  int ref_gnss_ch;
  const int visible_ch_num = common_observed_gnss_sat_id.size();
  Eigen::VectorXd h_sdcp = Eigen::VectorXd::Zero(visible_ch_num);
  Eigen::MatrixXd R_sdcp = R_.bottomRightCorner(visible_ch_num, visible_ch_num);
  Eigen::VectorXd Rv_sdcp(visible_ch_num);
  // 低仰角時におけるPCVの影響を小さくするためにelevationに閾値を設けることにする．
  const double ele_thresh_deg = 30; // deg
  int count = 0;
  std::vector<double> sdcp_raw;
  std::vector<int> gnss_ids;
  for (int ch = 0; ch < visible_ch_num; ch++)
  {
    const int gnss_id = common_observed_gnss_sat_id.at(ch);
    const int main_ch = GetIndexOfStdVector<int>(main_observation.info_.now_observed_gnss_sat_id, gnss_id);
    const int target_ch = GetIndexOfStdVector<int>(gnss_observations_.at(1).info_.now_observed_gnss_sat_id, gnss_id);

    double main_elevation_deg = main_observation.GetGnssElevationDeg(main_ch);
    if (main_elevation_deg < ele_thresh_deg) continue;

    const double pcc_main = x_est_main.pcc->GetPCC_m(main_observation.GetGnssAzimuthDeg(main_ch), main_elevation_deg);
    const double carrier_phase_main = main_observation.CalculateCarrierPhase(gnss_id, ConvEigenVecToLibraVec<3>(x_est_main.position), x_est_main.clock(0), x_est_main.ambiguity.N.at(main_ch), L1_lambda, pcc_main);

    const double pcc_target = x_est_target.pcc->GetPCC_m(gnss_observations_.at(1).GetGnssAzimuthDeg(target_ch), gnss_observations_.at(1).GetGnssElevationDeg(target_ch));
    const double carrier_phase_target = gnss_observations_.at(1).CalculateCarrierPhase(gnss_id, ConvEigenVecToLibraVec<3>(x_est_target.position), x_est_target.clock(0), x_est_target.ambiguity.N.at(target_ch), L1_lambda, pcc_target);

    h_sdcp(count) = carrier_phase_target - carrier_phase_main;
    sdcp_raw.push_back(sdcp_vec.at(ch));
    gnss_ids.push_back(gnss_id);
    Rv_sdcp(count) = R_sdcp(ch, ch);

  // 視線方向が重要であることを考えて仰角最大の衛星を参照にする．
    if (main_elevation_deg > max_elevation)
    {
      max_elevation = main_elevation_deg;
      ref_gnss_ch = count;
    }
    count++;
  }

  if (count == 0) return;

  h_sdcp.conservativeResize(count);
  Rv_sdcp.conservativeResize(count);
  R_sdcp = Rv_sdcp.asDiagonal();
  Eigen::MatrixXd M_dd = Eigen::MatrixXd::Zero(count - 1, count);
  Eigen::VectorXd res_ddcp = Eigen::VectorXd::Zero(count - 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(count - 1, 3);
  int ddcp_ch_offset = 0;
  for (int ch = 0; ch < count; ch++)
  {
    if (ch == ref_gnss_ch)
    {
      ddcp_ch_offset = -1;
      continue;
    }

    const double z_ddcp = sdcp_raw.at(ch) - sdcp_raw.at(ref_gnss_ch);
    const double h_ddcp = h_sdcp(ch) - h_sdcp(ref_gnss_ch);

    // std::cout << "z_ddcp" << z_ddcp << std::endl;
    // std::cout << "h_ddcp" << h_ddcp << std::endl;

    res_ddcp(ch + ddcp_ch_offset) = z_ddcp - h_ddcp;
    M_dd(ch + ddcp_ch_offset, ch) = 1.0; M_dd(ch + ddcp_ch_offset, ref_gnss_ch) = -1.0;
    // 視線方向ベクトルは2衛星でほぼ同等とみなしてmainのものを考える．
    const int main_ch = GetIndexOfStdVector(main_observation.info_.now_observed_gnss_sat_id, gnss_ids.at(ch));
    const int main_ref_gnss_original_ch = GetIndexOfStdVector(main_observation.info_.now_observed_gnss_sat_id, gnss_ids.at(ref_gnss_ch));
    const libra::Vector<3> de = main_observation.GetGnssDirection_c(main_ch) - main_observation.GetGnssDirection_c(main_ref_gnss_original_ch);

    for (int i = 0; i < 3; i++) H(ch + ddcp_ch_offset, i) = de[i];
    // H.block(ch, 0, 1, 3) = ConvLibraVecToEigenVec(h);
  }
  Eigen::MatrixXd R_ddcp = M_dd * R_sdcp * M_dd.transpose();

  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(R_ddcp.rows(), R_ddcp.cols());
  for (int i = 0; i < count - 1; i++)
  {
    // if (fabs(R_ddcp(i, i)) > 1e-10) W(i, i) = 1.0 / R_ddcp(i, i);
    for (int j = i; j < count - 1; j++)
    {
      // 0割りしないようにする．
      if (fabs(R_ddcp(i, j)) > 1e-10){ W(i, j) = 1.0 / R_ddcp(i, j); W(j, i) = W(i, j); }
      // if (i != j){ W(i, j) = 0.0; W(j, i) = 0.0; }
    }
  }
  // std::cout << "W" << W << std::endl;

  // W = Eigen::MatrixXd::Identity(R_ddcp.rows(), R_ddcp.cols());
  x_est_target.pcc->DpcoInitialEstimation(H, res_ddcp, W);
}

void PBD_dgps::GetStateFromVector(const int num_main_state_all, const Eigen::VectorXd& x_state)
{
  x_est_main.position = x_state.topRows(3);
  x_est_main.clock = x_state.block(3, 0, 1, 1);
  x_est_main.velocity = x_state.block(4, 0, 3, 1);

  x_est_target.position = x_state.block(num_main_state_all, 0, 3, 1);
  x_est_target.clock = x_state.block(num_main_state_all + 3, 0, 1, 1);
  x_est_target.velocity = x_state.block(num_main_state_all + 4, 0, 3, 1);
#ifdef REDUCED_DYNAMIC
  x_est_main.acceleration = x_state.block(7, 0, 3, 1);
  x_est_target.acceleration = x_state.block(num_main_state_all + 7, 0, 3, 1);
#endif // REDUCED_DYNAMIC
}
