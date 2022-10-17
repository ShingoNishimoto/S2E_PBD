#include <iomanip>
#include <set>
#include "PBD_dgps.h"
#include "PBD_Lambda.h"
#include <GeoPotential.h>

// clock model
#define CLOCK_IS_WHITE_NOISE (1)
// #define CLOCK_IS_RANDOM_WALK (1)

// inverse calculation method
#define CHOLESKY (0)
#define QR (1)

// Kalman Filter method
// #define AKF
// #define AIR_DRAG_ON

// #define REDUCED_DYNAMIC

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

#define N_DEBUG

#undef cross

static const int conv_index_from_gnss_sat_id(std::vector<int> observed_gnss_sat_id, const int gnss_sat_id);
static const double PBD_DGPS_kConvNm2m = 1e-9;
static const int precision = 7; // position

// outputを変えるときは"result.csv"を変更する．せめてパスは変えたい．
PBD_dgps::PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Dynamics& main_dynamics, const Dynamics& target_dynamics, PBD_GnssObservation& main_observation, PBD_GnssObservation& target_observation, PBD_GeoPotential* geop) :mt(42), step_time(sim_time_.GetStepSec()), ofs("result_new.csv"), num_of_gnss_satellites_(gnss_satellites_.GetNumOfSatellites()), main_dynamics_(main_dynamics), target_dynamics_(target_dynamics), receiver_clock_bias_main_(main_observation.receiver_clock_bias_), receiver_clock_bias_target_(target_observation.receiver_clock_bias_), gnss_observations_({ main_observation, target_observation }), geo_potential_(geop)
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

  x_est_.push_back(x_est_main);
  x_est_.push_back(x_est_target);
  // std::vector<PBD_GnssObservation&> gnss_observations_{ main_observation , target_observation }
  // ここは関係ない？
  // gnss_observations_.push_back(main_observation);
  // gnss_observations_.push_back(target_observation);

  GnssObserveModel main_model{};
  GnssObserveModel target_model{};
  gnss_observed_models_.push_back(main_model);
  gnss_observed_models_.push_back(target_model);

  antenna_pos_b_.push_back(main_observation.GetAntennaPosition());
  antenna_pos_b_.push_back(target_observation.GetAntennaPosition());

  true_N_main = Eigen::VectorXd::Zero(NUM_GNSS_CH);
  true_N_target = Eigen::VectorXd::Zero(NUM_GNSS_CH);

  // STMの初期化
  Phi_ = Eigen::Matrix<double, NUM_STATE_ALL, NUM_STATE_ALL>::Identity();
  Phi_(3, 3) = 0.0; Phi_(NUM_SINGLE_STATE_ALL + 3, NUM_SINGLE_STATE_ALL + 3) = 0.0;

  // 初期分散
  std::normal_distribution<> position_dist(0.0,sigma_r_ini);
  std::normal_distribution<> receiver_clock_dist(0.0, sigma_cdt_ini);
  std::normal_distribution<> velocity_dist(0.0, sigma_v_ini);
  std::normal_distribution<> acc_r_dist(0.0, sigma_acc_r_ini);
  std::normal_distribution<> acc_t_dist(0.0, sigma_acc_t_ini);
  std::normal_distribution<> acc_n_dist(0.0, sigma_acc_n_ini);
  std::normal_distribution<> N_dist(0.0, sigma_N_ini);

  // Vって名前は変えた方がいいかも
  Eigen::VectorXd V = Eigen::VectorXd::Constant(NUM_STATE_ALL, 0);

  // 状態量を減らした部分を実装．???????????????????????????????????????
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

  // 初期は全部に入れている．
  // for (int i = 0; i < NUM_GNSS_CH; ++i) V(NUM_SINGLE_STATE + i) = pow(sigma_N_ini, 2.0); // N

  // 以下がtarget
  for (int i = 0; i < 3; ++i) V(NUM_SINGLE_STATE_ALL + i) = pow(sigma_r_ini, 2.0);
  V(NUM_SINGLE_STATE_ALL + 3) = pow(sigma_cdt_ini, 2.0);
  for (int i = 0; i < 3; ++i) V(NUM_SINGLE_STATE_ALL + 4 + i) = pow(sigma_v_ini, 2.0);
#ifdef REDUCED_DYNAMIC
  V(NUM_SINGLE_STATE_ALL + 7) = pow(sigma_acc_r_ini, 2.0);
  V(NUM_SINGLE_STATE_ALL + 8) = pow(sigma_acc_t_ini, 2.0);
  V(NUM_SINGLE_STATE_ALL + 9) = pow(sigma_acc_n_ini, 2.0);
#endif // REDUCED_DYNAMIC

  // for (int i = 0; i < NUM_GNSS_CH; ++i) V(NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + i) = pow(sigma_N_ini, 2.0); // N

  M = V.asDiagonal(); // 誤差分散共分散行列M

  // Process noise Q
  CalculateQ();
  // Measurement noise R
  Eigen::VectorXd Rv = Eigen::VectorXd::Zero(NUM_OBSERVABLES);
  Rv.topRows(2*NUM_GNSS_CH) = (pow(pseudo_sigma / 2.0, 2.0) + pow(carrier_sigma / 2.0, 2.0))* Eigen::VectorXd::Ones(2*NUM_GNSS_CH); // GRAPHIC
  Rv.bottomRows(NUM_GNSS_CH) = 2.0*pow(carrier_sigma, 2.0) * Eigen::VectorXd::Ones(NUM_GNSS_CH); // SDCP
  R = Rv.asDiagonal();

  // 初期位置はガウシアンからサンプリング．mtは乱数のシード
  for(int i = 0; i < 3; ++i) x_est_main.position(i) += position_dist(mt);
  x_est_main.clock(0) += receiver_clock_dist(mt);
  for (int i = 0; i < 3; ++i) x_est_main.velocity(i) += velocity_dist(mt);

#ifdef REDUCED_DYNAMIC
  // x_est_main.acceleration(0) += acc_r_dist(mt);
  // x_est_main.acceleration(1) += acc_t_dist(mt);
  // x_est_main.acceleration(2) += acc_n_dist(mt);
#endif
  // for(int i = 0; i < NUM_GNSS_CH; ++i) x_est_main.ambiguity.N.at(i) += N_dist(mt);

  for(int i = 0; i < 3; ++i) x_est_target.position(i) += position_dist(mt);
  x_est_target.clock(0) += receiver_clock_dist(mt);
  for (int i = 0; i < 3; ++i) x_est_target.velocity(i) += velocity_dist(mt);

#ifdef REDUCED_DYNAMIC
  // x_est_target.acceleration(0) += acc_r_dist(mt);
  // x_est_target.acceleration(1) += acc_t_dist(mt);
  // x_est_target.acceleration(2) += acc_n_dist(mt);
#endif
  // for(int i = 0; i < NUM_GNSS_CH; ++i) x_est_target.ambiguity.N.at(i) += N_dist(mt);

  common_observed_status.assign(num_of_gnss_satellites_, false);
  // for (int i = 0; i < NUM_GNSS_CH; ++i) main_free_ch.push_back(i);
  // for (int i = 0; i < NUM_GNSS_CH; ++i) common_free_ch.push_back(i);

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
  // 参照渡ししたものを代入するのは無理．
  gnss_observations_.clear();
  gnss_observations_.push_back(main_observation_);
  gnss_observations_.push_back(target_observation_);

  // これがずれている気がするな．
  trans_eci_to_ecef_ = earth_rotation.GetDCMJ2000toXCXF();

  double elapsed_time = sim_time_.GetElapsedSec();
  double tmp = floor(elapsed_time/observe_step_time + 1e-4); //1e-4は数値誤差
  double tmp_log = floor(elapsed_time/log_step_time + 1e-4);

  //まず更新
  OrbitPropagation();

  //観測時間にピッタリ
  if(abs(elapsed_time - tmp*observe_step_time) < step_time/2.0){
    // M = UpdateM(); // これを入れるのがどこであるべきなのか?

    // ここの引数をどう渡すかとかが関係している？
    // PBD_GnssObservation& main_observation_ = gnss_observations_.at(0);
    SetBiasToObservation(0, x_est_main, gnss_observations_.at(0));
    // PBD_GnssObservation& target_observation_ = gnss_observations_.at(1);
    SetBiasToObservation(1, x_est_target, gnss_observations_.at(1));

    KalmanFilter();
  }

  //log output
  if(abs(elapsed_time - tmp_log*log_step_time) < step_time/2.0){
    libra::Vector<3> sat_position_main = main_dynamics_.GetPosition_i();
    libra::Vector<3> sat_velocity_main = main_dynamics_.GetOrbit().GetSatVelocity_i();
    libra::Vector<3> sat_position_target = target_dynamics_.GetPosition_i();
    libra::Vector<3> sat_velocity_target = target_dynamics_.GetOrbit().GetSatVelocity_i();

    // ECIでの真値（位置，クロックバイアス，速度）を残す．
    for(int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << sat_position_main[i] << ","; // r_m_true
    ofs << std::fixed << std::setprecision(precision) << receiver_clock_bias_main_ << ","; // t_m_true
    for(int i = 0;i < 3;++i) ofs << std::fixed << std::setprecision(precision) << sat_velocity_main[i] << ","; // v_m_true

    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << sat_position_target[i] << ","; // r_t_true
    ofs << std::fixed << std::setprecision(precision) << receiver_clock_bias_target_ << ","; // t_t_true
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
    Eigen::Matrix3d trans_eci_to_rtn_main = TransRTN2ECI(sat_pos_eci_main, sat_vel_eci_main).inverse();
    Eigen::Matrix3d trans_eci_to_rtn_target = TransRTN2ECI(sat_pos_eci_target, sat_vel_eci_target).inverse();

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

    for (int i = 0; i < NUM_GNSS_CH; ++i) ofs << std::fixed << std::setprecision(precision) << true_N_main(i) << ","; // N_true
    for (int i = 0; i < NUM_GNSS_CH; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_main.ambiguity.N.at(i) << ","; // N_est
    for (int i = 0; i < NUM_GNSS_CH; ++i) ofs << std::fixed << std::setprecision(precision) << true_N_target(i) << ","; // N_true
    for (int i = 0; i < NUM_GNSS_CH; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_target.ambiguity.N.at(i) << ","; // N_est
    // FIXME: ここもRTNに変換する．
    for (int i = 0; i < NUM_STATE_ALL; ++i) ofs << std::fixed << std::setprecision(precision) << M(i, i) << ",";
    // record visible gnss sat number
    /*
    int ans = 0;
    for(int i = 0;i < num_of_gnss_satellites_;++i){
        auto gnss_position = gnss_satellites_.GetSatellitePositionEci(i);
        if(CheckCanSeeSatellite(sat_position, gnss_position)) ++ans;
    }
    ofs << ans << std::endl;
    */
    // そもそもここでログをとるのが適切ではない．
    int visible_gnss_num_main = gnss_observations_.at(0).info_.now_observed_gnss_sat_id.size();
    int visible_gnss_num_target = gnss_observations_.at(1).info_.now_observed_gnss_sat_id.size();
    int visible_gnss_num_common = common_observed_gnss_sat_id.size();
    ofs << visible_gnss_num_main << ",";
    ofs << visible_gnss_num_target << ",";
    ofs << visible_gnss_num_common << ",";
    // main observe gnss sat id
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i >= visible_gnss_num_main) ofs << -1 << ",";
      else ofs << gnss_observations_.at(0).info_.now_observed_gnss_sat_id.at(i) << ",";
    }
    // target observe gnss sat id
    for (int i = 0; i < NUM_GNSS_CH; ++i)
    {
      if (i >= visible_gnss_num_target) ofs << -1 << ",";
      else ofs << gnss_observations_.at(1).info_.now_observed_gnss_sat_id.at(i) << ",";
    }
    // Q
    for (int i = 0; i < NUM_STATE_ALL; i++) ofs << Q(i, i) << ","; // 全部残す意味はないのかも
    // R
    for (int i = 0; i < NUM_OBSERVABLES; i++) ofs << R(i, i) << ",";

    // acc eci
    Eigen::Vector3d acc_m_i = TransRTN2ECI(x_est_main.position, x_est_main.velocity) * x_est_main.acceleration; // [nm/s2]
    Eigen::Vector3d acc_t_i = TransRTN2ECI(x_est_target.position, x_est_target.velocity) * x_est_target.acceleration; // [nm/s2]
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << acc_m_i(i) << ","; // a_m_i
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << acc_t_i(i) << ","; // a_t_i
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_main.acc_dist(i) << ","; // a_disturbance_m
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(precision) << x_est_target.acc_dist(i) << ","; // a_disturbance_t
    ofs << std::endl;
  }

  return;
}

void PBD_dgps::InitLogTable(void)
{

}

void PBD_dgps::OrbitPropagation()
{
  //RK4
  Eigen::MatrixXd Phi_main = Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE>::Identity();
  Eigen::MatrixXd Phi_target = Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE>::Identity();
  Phi_main(3, 3) = 0.0; Phi_target(3, 3) = 0.0;

  RK4(x_est_main.position, x_est_main.velocity, x_est_main.acceleration, x_est_main.acc_dist, Phi_main);
  Phi_.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Phi_main;

  RK4(x_est_target.position, x_est_target.velocity, x_est_target.acceleration, x_est_target.acc_dist, Phi_target);
  Phi_.block(NUM_SINGLE_STATE_ALL, NUM_SINGLE_STATE_ALL, NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Phi_target;

#ifdef REDUCED_DYNAMIC
  // acceleration ここでは q = sigma_acc_process
  Eigen::MatrixXd Phi_a = CalculatePhi_a(step_time); // ここもチューニングできないと意味ないのでは？
  double phi = Phi_a(0, 0);
  // double phi_t = Phi_a(1, 1);
  // double phi_n = Phi_a(2, 2);
  x_est_main.acceleration = Phi_a * x_est_main.acceleration;
  x_est_target.acceleration = Phi_a * x_est_target.acceleration;
  std::normal_distribution<> acc_r_process_noise(0.0, sigma_acc_r_process * sqrt(1 - phi * phi));
  std::normal_distribution<> acc_t_process_noise(0.0, sigma_acc_t_process * sqrt(1 - phi * phi));
  std::normal_distribution<> acc_n_process_noise(0.0, sigma_acc_n_process*sqrt(1 - phi*phi));
  // for (int i = 0; i < 3; ++i) x_est_main.acceleration(i) += acc_r_process_noise(mt);
  // for (int i = 0; i < 3; ++i) x_est_target.acceleration(i) += acc_r_process_noise(mt);
  x_est_main.acceleration(0) += acc_r_process_noise(mt);
  x_est_target.acceleration(0) += acc_r_process_noise(mt);
  x_est_main.acceleration(1) += acc_t_process_noise(mt);
  x_est_target.acceleration(1) += acc_t_process_noise(mt);
  x_est_main.acceleration(2) += acc_n_process_noise(mt);
  x_est_target.acceleration(2) += acc_n_process_noise(mt);
#endif // REDUCED_DYNAMIC

  // cdt
#ifdef CLOCK_IS_RANDOM_WALK
  std::normal_distribution<> cdt_process_noise(0.0, sigma_cdt_process*sqrt(step_time/tau_cdt));
  std::normal_distribution<> cdt_process_noise(0.0, sigma_cdt_process);
  x_est_main.clock(0) = cdt_process_noise(mt);
  x_est_target.clock(0) = cdt_process_noise(mt);
#endif // CLOCK_IS_RANDOM_WALK
  M = UpdateM();
  // M = Phi_*M*Phi_.transpose(); // プロセスノイズなし
}

void PBD_dgps::RK4(Eigen::Vector3d& position, Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration, Eigen::Vector3d& acc_dist, Eigen::MatrixXd& Phi)
{
  Eigen::Vector3d k0 = PositionDifferential(velocity);
  Eigen::Vector3d l0 = VelocityDifferential(position, velocity, acceleration, acc_dist);
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> n0 = CalculateJacobian(position, velocity);

  Eigen::Vector3d tmp_position = position + k0 * step_time / 2.0;
  Eigen::Vector3d tmp_velocity = velocity + l0 * step_time / 2.0;
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> tmp_Phi = Phi + n0 * step_time / 2.0;
  Eigen::Vector3d k1 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l1 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> n1 = CalculateJacobian(tmp_position, tmp_velocity);

  tmp_position = position + k1 * step_time / 2.0;
  tmp_velocity = velocity + l1 * step_time / 2.0;
  tmp_Phi = Phi + n1 * step_time / 2.0;
  Eigen::Vector3d k2 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l2 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> n2 = CalculateJacobian(tmp_position, tmp_velocity);

  tmp_position = position + k2 * step_time;
  tmp_velocity = velocity + l2 * step_time;
  tmp_Phi = Phi + n2 * step_time;
  Eigen::Vector3d k3 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l3 = VelocityDifferential(tmp_position, tmp_velocity, acceleration, acc_dist);
  Eigen::Matrix<double, NUM_SINGLE_STATE, NUM_SINGLE_STATE> n3 = CalculateJacobian(tmp_position, tmp_velocity);

  position += step_time * (k0 + 2.0 * k1 + 2.0 * k2 + k3) / 6.0;
  velocity += step_time * (l0 + 2.0 * l1 + 2.0 * l2 + l3) / 6.0;
  Phi      += step_time * (n0 + 2.0 * n1 + 2.0 * n2 + n3) / 6.0;
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

  // AddGeoPotentialDisturbance(position, acc_dist);

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

Eigen::MatrixXd PBD_dgps::UpdateM()
{
  // clock
#ifdef CLOCK_IS_RANDOM_WALK
  Phi_(3, 3) = 1.0; // random walkなのでΦは1
  Phi_(NUM_SINGLE_STATE_ALL + 3, NUM_SINGLE_STATE_ALL + 3) = 1.0;
#endif // CLOCK_IS_RANDOM_WALK
#ifdef CLOCK_IS_WHITE_NOISE
  Phi_(3, 3) = 0.0;
  Phi_(NUM_SINGLE_STATE_ALL + 3, NUM_SINGLE_STATE_ALL + 3) = 0.0;
#endif // CLOCK_IS_WHITE_NOISE

  // a
#ifdef REDUCED_DYNAMIC
  Phi_.block(7, 7, 3, 3) = CalculatePhi_a(observe_step_time);
  Phi_.block(NUM_SINGLE_STATE_ALL + 7, NUM_SINGLE_STATE_ALL + 7, 3, 3) = CalculatePhi_a(observe_step_time);
#endif // REDUCED_DYNAMIC

// N
Phi_.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, NUM_GNSS_CH, NUM_GNSS_CH) = Eigen::MatrixXd::Identity(NUM_GNSS_CH, NUM_GNSS_CH);
Phi_.bottomRightCorner(NUM_GNSS_CH, NUM_GNSS_CH) = Eigen::MatrixXd::Identity(NUM_GNSS_CH, NUM_GNSS_CH);

#ifdef AKF
  // ここでは更新しない．更新しないつもりなのに初期値に依存しているってことは更新されてしまっている？
#else
  CalculateQ();
#endif // AKF

  int n_main = gnss_observations_.at(0).info_.now_observed_gnss_sat_id.size();
  Q.block(NUM_SINGLE_STATE + n_main, NUM_SINGLE_STATE + n_main, NUM_GNSS_CH - n_main, NUM_GNSS_CH - n_main) = Eigen::MatrixXd::Zero(NUM_GNSS_CH - n_main, NUM_GNSS_CH - n_main);
  int n_target = gnss_observations_.at(1).info_.now_observed_gnss_sat_id.size();
  Q.block(NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + n_target, NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + n_target, NUM_GNSS_CH - n_target, NUM_GNSS_CH - n_target) = Eigen::MatrixXd::Zero(NUM_GNSS_CH - n_target, NUM_GNSS_CH - n_target);

  Eigen::MatrixXd res = Phi_ * M * Phi_.transpose() + Q;
  // Nのない部分を0に落とす <- これが必要なのかもあやしい．
  res.block(NUM_SINGLE_STATE + n_main, NUM_SINGLE_STATE + n_main, NUM_GNSS_CH - n_main, NUM_GNSS_CH - n_main) = Eigen::MatrixXd::Zero(NUM_GNSS_CH - n_main, NUM_GNSS_CH - n_main);

  res.block(NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + n_target, NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + n_target, NUM_GNSS_CH - n_target, NUM_GNSS_CH - n_target) = Eigen::MatrixXd::Zero(NUM_GNSS_CH - n_target, NUM_GNSS_CH - n_target);

  return res;
}

Eigen::MatrixXd PBD_dgps::CalculateA(const Eigen::Vector3d& position_main, const Eigen::Vector3d& velocity_main, const Eigen::Vector3d& position_target, const Eigen::Vector3d& velocity_target)
{
  // Dynamics and Kinematics model
  Eigen::MatrixXd Jacobi_main = CalculateJacobian(position_main, velocity_main);
  Eigen::MatrixXd Jacobi_target = CalculateJacobian(position_target, velocity_target);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(NUM_STATE_ALL, NUM_STATE_ALL);
  A.topLeftCorner(NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Jacobi_main;
  A.block(NUM_SINGLE_STATE_ALL, NUM_SINGLE_STATE_ALL, NUM_SINGLE_STATE, NUM_SINGLE_STATE) = Jacobi_target;

  return A;
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
  // J2
  // double J2_coefficient = 3.0 / 2.0 * mu_e * J2_const * pow(Earth_Radius, 2.0);
  // libra::Matrix<3, 3> A_j2_ecef(0);
  // A_j2_ecef[0][0] = - J2_coefficient / pow(r, 5.0) * (1.0 - 5.0 * (pow(x_ef, 2.0) + pow(z_ef, 2.0)) / pow(r, 2.0) + 35.0 * pow(x_ef * z_ef, 2.0) / pow(r, 4.0));
  // A_j2_ecef[0][1] = - J2_coefficient * x_ef * y_ef / pow(r, 7.0) * (-5.0 + 35.0 * pow(z_ef / r, 2.0));
  // A_j2_ecef[0][2] = - J2_coefficient * x_ef * z_ef / pow(r, 7.0) * (-15.0 + 35.0 * pow(z_ef / r, 2.0));
  // A_j2_ecef[1][0] = - J2_coefficient * x_ef * y_ef / pow(r, 7.0) * (-5.0 + 35.0 * pow(z_ef / r, 2.0));
  // A_j2_ecef[1][1] = - J2_coefficient / pow(r, 5.0) * (1.0 - 5.0 * (pow(y_ef, 2.0) + pow(z_ef, 2.0)) / pow(r, 2.0) + 35.0 * pow(y_ef * z_ef, 2.0) / pow(r, 4.0));
  // A_j2_ecef[1][2] = - J2_coefficient * y_ef * z_ef / pow(r, 7.0) * (-15.0 + 35.0 * pow(z_ef / r, 2.0));
  // A_j2_ecef[2][0] = - J2_coefficient * x_ef * z_ef / pow(r, 7.0) * (-15.0 + 35.0 * pow(z_ef / r, 2.0));
  // A_j2_ecef[2][1] = - J2_coefficient * y_ef * z_ef / pow(r, 7.0) * (-15.0 + 35.0 * pow(z_ef / r, 2.0));
  // A_j2_ecef[2][2] = - J2_coefficient / pow(r, 5.0) * (3.0 - 30.0 * pow(z_ef/r, 2.0) + 35.0 * pow(z_ef/r, 4.0));

  // libra::Matrix<3, 3> A_j2_eci = libra::transpose(trans_eci_to_ecef_) * A_j2_ecef * trans_eci_to_ecef_;

  // for (uint8_t i = 0; i < 3; i++)
  // {
  //   for (uint8_t j = 0; j < 3; j++) A(i + 4, j) += A_j2_eci[i][j];
  // }

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
  Eigen::Vector3d rtn_r = position.normalized();
  Eigen::Vector3d rtn_n = position.cross(velocity);
  rtn_n.normalize();
  Eigen::Vector3d rtn_t = rtn_n.cross(rtn_r);
  rtn_t.normalize();

  Eigen::MatrixXd RTN2ECI(3,3);
  RTN2ECI.block(0, 0, 3, 1) = rtn_r;
  RTN2ECI.block(0, 1, 3, 1) = rtn_t;
  RTN2ECI.block(0, 2, 3, 1) = rtn_n;
  return RTN2ECI;
};

// tに対するノイズの加え方が間違っている気がする．これのせいか？
void PBD_dgps::CalculateQ(void)
{
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(NUM_STATE_ALL, 8); // aとclock
  B(3, 0) = 1.0;
  B(NUM_SINGLE_STATE_ALL + 3, 4) = 1.0;
#ifdef REDUCED_DYNAMIC
  B.block(7, 1, 3, 3) = Eigen::Matrix3d::Identity();
  B.block(NUM_SINGLE_STATE_ALL + 7, 5, 3, 3) = Eigen::Matrix3d::Identity();
#endif // REDUCED_DYNAMIC

  Eigen::MatrixXd Q_at = CalculateQ_at();
  // 観測するたびにNの部分を初期化？

  // Q = BQ_atB^t
  Q = B * Q_at * B.transpose();
  // add process noise for r and v
  Q.block(0, 0, 3, 3) = pow(sigma_r_process, 2.0) * Eigen::Matrix3d::Identity() * pow(step_time, 2.0);
  Q.block(4, 4, 3, 3) = pow(sigma_v_process, 2.0) * Eigen::Matrix3d::Identity() * pow(step_time, 2.0);
  Q.block(NUM_SINGLE_STATE_ALL, NUM_SINGLE_STATE_ALL, 3, 3) = pow(sigma_r_process, 2.0) * Eigen::Matrix3d::Identity() * pow(step_time, 2.0);
  Q.block(NUM_SINGLE_STATE_ALL + 4, NUM_SINGLE_STATE_ALL + 4, 3, 3) = pow(sigma_v_process, 2.0) * Eigen::Matrix3d::Identity() * pow(step_time, 2.0);

#ifndef N_DEBUG
  // N process
  Q.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, NUM_GNSS_CH, NUM_GNSS_CH) = pow(sigma_N_process, 2.0) * Eigen::MatrixXd::Identity(NUM_GNSS_CH, NUM_GNSS_CH) *pow(step_time, 2.0);
  Q.block(NUM_SINGLE_STATE + NUM_SINGLE_STATE_ALL, NUM_SINGLE_STATE + NUM_SINGLE_STATE_ALL, NUM_GNSS_CH, NUM_GNSS_CH) = pow(sigma_N_process, 2.0) * Eigen::MatrixXd::Identity(NUM_GNSS_CH, NUM_GNSS_CH) * pow(step_time, 2.0);
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
  double phi = exp(-step_time / tau_a); // ここはobserved_timeとどっちなのか？
  q_acc_r = pow(sigma_acc_r_process, 2.0) * (1 - pow(phi, 2.0));
  q_acc_t = pow(sigma_acc_t_process, 2.0) * (1 - pow(phi, 2.0));
  q_acc_n = pow(sigma_acc_n_process, 2.0) * (1 - pow(phi, 2.0));
#endif // REDUCED_DYNAMIC

  double q_cdt;
#ifdef CLOCK_IS_WHITE_NOISE
  q_cdt = pow(sigma_cdt_process * step_time, 2.0);
#endif // CLOCK_IS_WHITE_NOISE
#ifdef CLOCK_IS_RANDOM_WALK
  q_cdt = pow(sigma_cdt_process, 2.0) * (step_time / tau_cdt);
#endif // CLOCK_IS_RANDOM_WALK

  Q_at(0, 0) = q_cdt;
  Q_at(1, 1) = q_acc_r; Q_at(2, 2) = q_acc_t; Q_at(3, 3) = q_acc_n;

  Q_at(4, 4) = q_cdt;
  Q_at(5, 5) = q_acc_r; Q_at(6, 6) = q_acc_t; Q_at(7, 7) = q_acc_n;

  return Q_at;
}

Eigen::MatrixXd PBD_dgps::CalculatePhi_a(const double dt)
{
  // ここが各軸に応じて変わる部分なのか．
  double phi = exp(-dt / tau_a); // constやないか
  Eigen::MatrixXd Phi_a = Eigen::MatrixXd::Identity(3, 3);
#ifdef REDUCED_DYNAMIC
  Phi_a *= phi;
#endif // REDUCED_DYNAMIC

  return Phi_a;
};

void PBD_dgps::KalmanFilter()
{
  // gnss_observed_hoge の観測される時間はどうなっているのか？

  // GRAPHIC*2 + SDCPにする．
  Eigen::VectorXd z = Eigen::VectorXd::Zero(NUM_OBSERVABLES); //観測ベクトル
  Eigen::VectorXd h_x = Eigen::VectorXd::Zero(NUM_OBSERVABLES); // 観測モデル行列

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(NUM_OBSERVABLES, NUM_STATE_ALL); //観測行列（dhx/dx）
  Eigen::VectorXd R_V = Eigen::MatrixXd::Zero(NUM_OBSERVABLES, 1); // 観測誤差共分散．
  UpdateObservations(z, h_x, H, R_V);

  int n_main = gnss_observations_.at(0).info_.now_observed_gnss_sat_id.size();
  int offset = n_main;
  int size = NUM_GNSS_CH - n_main;
  int n_target = gnss_observations_.at(1).info_.now_observed_gnss_sat_id.size();
  int n_common = common_observed_gnss_sat_id.size();

// #ifdef AKF
//   // Rは更新しない．
//   R.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);
//   offset = NUM_GNSS_CH + n_target;
//   size = NUM_GNSS_CH - n_target;
//   R.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);
//   offset = 2*NUM_GNSS_CH + n_common;
//   size = NUM_GNSS_CH - n_common;
//   R.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);
// #else
  R = R_V.asDiagonal();
// #endif

  Eigen::MatrixXd hmh = H * M * H.transpose();

  Eigen::MatrixXd tmp = R + hmh; // (observation_num, observation_num)
  // 観測量がないところは0にし直すみたいな加工が必要かもしれない．数字の桁数で打ち切りみたいな形にできればいい．
  //カルマンゲイン
  Eigen::MatrixXd K = CalculateK(H, tmp);

  Eigen::VectorXd x_predict = Eigen::VectorXd(NUM_STATE_ALL);
  x_predict.topRows(3) = x_est_main.position;
  x_predict.block(3, 0, 1, 1) = x_est_main.clock;
  x_predict.block(4, 0, 3, 1)  = x_est_main.velocity;
#ifdef REDUCED_DYNAMIC
  x_predict.block(7, 0, 3, 1) = x_est_main.acceleration;
  x_predict.block(NUM_SINGLE_STATE_ALL + 7, 0, 3, 1) = x_est_target.acceleration;
#endif // REDUCED_DYNAMIC
  x_predict.block(NUM_SINGLE_STATE, 0, NUM_GNSS_CH, 1) = Eigen::Map<Eigen::VectorXd>(&x_est_main.ambiguity.N[0], x_est_main.ambiguity.N.size());
  //x_predict(10) = Cd;
  x_predict.block(NUM_SINGLE_STATE_ALL, 0, 3, 1) = x_est_target.position;
  x_predict.block(NUM_SINGLE_STATE_ALL + 3, 0, 1, 1) = x_est_target.clock;
  x_predict.block(NUM_SINGLE_STATE_ALL + 4, 0, 3, 1) = x_est_target.velocity;
  x_predict.bottomRows(NUM_GNSS_CH) = Eigen::Map<Eigen::VectorXd>(&x_est_target.ambiguity.N[0], x_est_target.ambiguity.N.size());

  // innovationの記号を何にするかは要検討
  Eigen::VectorXd E_pre = z - h_x;
  // まずアンテナ位置に変換
  Eigen::VectorXd x_ant_predict = ConvReceivePosToCenterOfMass(x_predict);
  Eigen::VectorXd x_update = x_ant_predict + K * E_pre;
  // 重心位置に戻す．
  x_update = ConvReceivePosToCenterOfMass(x_update); // 参照渡しにしてもいいかも

 //更新
  x_est_main.position = x_update.topRows(3);
  x_est_main.clock = x_update.block(3, 0, 1, 1);
  x_est_main.velocity = x_update.block(4, 0, 3, 1);
  //Cd = x_update(10);
  x_est_target.position = x_update.block(NUM_SINGLE_STATE_ALL, 0, 3, 1);
  x_est_target.clock = x_update.block(NUM_SINGLE_STATE_ALL + 3, 0, 1, 1);
  x_est_target.velocity = x_update.block(NUM_SINGLE_STATE_ALL + 4, 0, 3, 1);
#ifdef REDUCED_DYNAMIC
  x_est_main.acceleration = x_update.block(7, 0, 3, 1);
  x_est_target.acceleration = x_update.block(NUM_SINGLE_STATE_ALL + 7, 0, 3, 1);
#endif // REDUCED_DYNAMIC
#ifndef N_DEBUG
  for (int i = 0; i < NUM_GNSS_CH; ++i)
  {
    x_est_main.ambiguity.N.at(i) = x_update(NUM_SINGLE_STATE + i);
    x_est_target.ambiguity.N.at(i) = x_update(NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + i);
  }
#else
    // For debug 一旦整数不定性は0とする．
    x_est_main.ambiguity.N = std::vector<double>(NUM_GNSS_CH, 0);
    x_est_target.ambiguity.N = std::vector<double>(NUM_GNSS_CH, 0);
#endif // N_DEBUG

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(NUM_STATE_ALL, NUM_STATE_ALL);
  tmp = (I - K*H);
  M = tmp*M*tmp.transpose() + K*R*K.transpose(); // 負になるのが存在しているが，これは何？

  // 数値誤差で発散してしまったやつを処理．
  offset = NUM_SINGLE_STATE + n_main;
  size = NUM_GNSS_CH - n_main;
  M.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);
  Q.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);
  offset = NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + n_target;
  size = NUM_GNSS_CH - n_target;
  M.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);
  Q.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);

#ifdef AKF
  // residual
  // GRAPHIC*2 + SDCPにする．
  Eigen::VectorXd h_x_post = Eigen::VectorXd::Zero(NUM_OBSERVABLES); // 観測モデル行列

  UpdateObservations(z, h_x_post, H, R_V);
  Eigen::VectorXd E_post = z - h_x_post;
  // FIXME: ここでSDCPの精度がめっちゃ悪いことになってしまっているのが原因な気がする．
  R = alpha * R + (1 - alpha) * (E_post*E_post.transpose() + H*M*H.transpose()); // residual based R-adaptation
  Eigen::MatrixXd Q_dash = alpha * Q + (1 - alpha) * K * E_pre * (K * E_pre).transpose(); // Innovation based Q-adaptation

  DynamicNoiseScaling(Q_dash, CalculateA(), H);
  // これしてるから，絶対軌道精度の影響（つまり残差の大きさ）ではなくて，収束してしまう？
/*
  Q = Eigen::MatrixXd::Zero(NUM_STATE_ALL, NUM_STATE_ALL);
  // a, cdtだけもらう．
  Q(3, 3) = Q_dash(3, 3);
  Q.block(7, 7, 3, 3) = Q_dash.block(7, 7, 3, 3);
  Q(NUM_SINGLE_STATE_ALL + 3, NUM_SINGLE_STATE_ALL + 3) = Q_dash(NUM_SINGLE_STATE_ALL + 3, NUM_SINGLE_STATE_ALL + 3);
  Q.block(NUM_SINGLE_STATE_ALL + 7, NUM_SINGLE_STATE_ALL + 7, 3, 3) = Q_dash.block(NUM_SINGLE_STATE_ALL + 7, NUM_SINGLE_STATE_ALL + 7, 3, 3);
*/
#endif // AKF

  // TODO: Nについて，収束している？ものは0に落とす．
  for (int i = 0; i < NUM_GNSS_CH; ++i)
  {
    int N_offset = NUM_SINGLE_STATE + i;
    if (M(N_offset, N_offset) < 0) M(N_offset, N_offset) = pow(sigma_N_ini, 2.0); // sqrt(pow(M(N_offset, N_offset), 2.0));
    N_offset = NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + i;
    if (M(N_offset, N_offset) < 0) M(N_offset, N_offset) = pow(sigma_N_ini, 2.0); // sqrt(pow(M(N_offset, N_offset), 2.0));
  }

  // others
  // for (int i = 0; i < NUM_SINGLE_STATE; ++i)
  // {
  //   int offset = i;
  //   if (M(offset, offset) < 1e-6) M(offset, offset) = sqrt(pow(M(offset, offset), 2.0));
  //   offset = NUM_SINGLE_STATE_ALL + i;
  //   if (M(offset, offset) < 1e-6) M(offset, offset) = sqrt(pow(M(offset, offset), 2.0));
  // }

  // TODO: make the Double Differential ambiguity

  // IAR

#ifdef LAMBDA_DEBUG
  Eigen::MatrixXd Q_a_main = M.block(NUM_SINGLE_STATE, NUM_SINGLE_STATE, NUM_GNSS_CH, NUM_GNSS_CH);
  Eigen::VectorXd a_main = x_est_main.N; // ambiguity
  // non-zeroに限定する．
  int num_of_obsetved_gnss =  gnss_observations_.at(0).info_.now_observed_gnss_sat_id.size();
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
    gnss_observed_models_.at(i).geometric_range.clear();
    gnss_observed_models_.at(i).pseudo_range_model.clear();
    gnss_observed_models_.at(i).carrier_phase_range_model.clear();
  }
  return;
}

Eigen::MatrixXd PBD_dgps::CalculateK(Eigen::MatrixXd H, Eigen::MatrixXd S)
{
  int observe_gnss_num_m = gnss_observations_.at(0).info_.now_observed_gnss_sat_id.size();
  int observe_gnss_num_t = gnss_observations_.at(1).info_.now_observed_gnss_sat_id.size();
  int observe_gnss_num_c = common_observed_gnss_sat_id.size();

  Eigen::MatrixXd MHt = M * H.transpose();
  ResizeS(S, observe_gnss_num_m, observe_gnss_num_t, observe_gnss_num_c);
  ResizeMHt(MHt, observe_gnss_num_m, observe_gnss_num_t, observe_gnss_num_c);
#if CHOLESKY
  Eigen::MatrixXd ST = S.transpose();
  Eigen::LDLT<Eigen::MatrixXd> LDLTOftmpT(ST);
  Eigen::MatrixXd KT = LDLTOftmpT.solve(MHt.transpose());
  Eigen::MatrixXd K_calc = KT.transpose();
#elif QR
  Eigen::MatrixXd ST = S.transpose();
  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QROftmpT(ST);
  Eigen::MatrixXd KT = QROftmpT.solve(MHt.transpose());
  Eigen::MatrixXd K_calc = KT.transpose();
#else
  Eigen::MatrixXd K_calc = MHt * S.inverse();
#endif

  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(NUM_STATE_ALL, NUM_OBSERVABLES);
  K.topLeftCorner(NUM_SINGLE_STATE + observe_gnss_num_m, observe_gnss_num_m) = K_calc.topLeftCorner(NUM_SINGLE_STATE + observe_gnss_num_m, observe_gnss_num_m);
  K.block(NUM_SINGLE_STATE_ALL, 0, NUM_SINGLE_STATE + observe_gnss_num_t, observe_gnss_num_m) = K_calc.block(NUM_SINGLE_STATE + observe_gnss_num_m, 0, NUM_SINGLE_STATE + observe_gnss_num_t, observe_gnss_num_m);
  K.block(0, NUM_GNSS_CH, NUM_SINGLE_STATE + observe_gnss_num_m, observe_gnss_num_t) = K_calc.block(0, observe_gnss_num_m, NUM_SINGLE_STATE + observe_gnss_num_m, observe_gnss_num_t);
  K.block(NUM_SINGLE_STATE_ALL, NUM_GNSS_CH, NUM_SINGLE_STATE + observe_gnss_num_t, observe_gnss_num_t) = K_calc.block(NUM_SINGLE_STATE + observe_gnss_num_m, observe_gnss_num_m, NUM_SINGLE_STATE + observe_gnss_num_t, observe_gnss_num_t);
  K.block(0, 2*NUM_GNSS_CH, NUM_SINGLE_STATE + observe_gnss_num_m, observe_gnss_num_c) = K_calc.block(0, observe_gnss_num_m+observe_gnss_num_t, NUM_SINGLE_STATE + observe_gnss_num_m, observe_gnss_num_c);
  K.block(NUM_SINGLE_STATE_ALL, 2*NUM_GNSS_CH, NUM_SINGLE_STATE + observe_gnss_num_t, observe_gnss_num_c) = K_calc.block(NUM_SINGLE_STATE + observe_gnss_num_m, observe_gnss_num_m+observe_gnss_num_t, NUM_SINGLE_STATE + observe_gnss_num_t, observe_gnss_num_c);
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

void PBD_dgps::UpdateObservationsGRAPHIC(const int sat_id, EstimatedVariables& x_est, const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
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
  double carrier_phase_range = (carrier_phase.first + carrier_phase.second) * L1_lambda;

  libra::Vector<3> sat_position{ 0.0 };
  for (int i = 0; i < 3; i++) sat_position[i] = x_est.position(i);
  const double sat_clock = x_est.clock(0);
  // ここら辺はGnssObserveModel内に格納する．
  double geometric_range = gnss_observation.CalculateGeometricRange(sat_position, gnss_position);
  observe_model.geometric_range.push_back(geometric_range);

  double pseudo_range_model = gnss_observation.CalculatePseudoRange(sat_position, gnss_position, sat_clock, gnss_clock);
  observe_model.pseudo_range_model.push_back(pseudo_range_model);

  double carrier_phase_range_model = gnss_observation.CalculateCarrierPhase(sat_position, gnss_position, sat_clock, gnss_clock, x_est.ambiguity.N.at(index), L1_lambda);
  observe_model.carrier_phase_range_model.push_back(carrier_phase_range_model);

  // GRAPHIC
  const int row_offset = sat_id*NUM_GNSS_CH + index;
  const int col_offset = sat_id*NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + index;
  z(row_offset) = (pseudo_range + carrier_phase_range) / 2;
  h_x(row_offset) = (pseudo_range_model + carrier_phase_range_model) / 2;
  for (int j = 0; j < 3; ++j) {
    const int col_index = sat_id * NUM_SINGLE_STATE_ALL + j;
    // position
    H(row_offset, col_index) = (x_est.position(j) - gnss_position[j]) / geometric_range;
  }
  // clock
  H(row_offset, sat_id * NUM_SINGLE_STATE_ALL + 3) = 1.0;
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
  const int row_offset = 2 * NUM_GNSS_CH + common_index;
  PBD_GnssObservation& main_observation = gnss_observations_.at(0);
  PBD_GnssObservation& target_observation = gnss_observations_.at(1);

  const int main_index = conv_index_from_gnss_sat_id(main_observation.info_.now_observed_gnss_sat_id, gnss_sat_id);
  // const int main_col_offset = 2 * NUM_GNSS_CH + common_index;
  const int target_index = conv_index_from_gnss_sat_id(target_observation.info_.now_observed_gnss_sat_id, gnss_sat_id);
  const int col_offset_main = NUM_SINGLE_STATE + main_index;
  const int col_offset_target = NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + target_index;

  // とりあえずL1を使う．
  // main
  const auto& carrier_phase_main = main_observation.observed_values_.L1_carrier_phase.at(main_index);
  double carrier_phase_range_main = (carrier_phase_main.first + carrier_phase_main.second) * x_est_main.lambda; // FIXME: x_est_をmaster情報にできるように修正する．
  // target
  const auto& carrier_phase_target = target_observation.observed_values_.L1_carrier_phase.at(target_index);
  double carrier_phase_range_target = (carrier_phase_target.first + carrier_phase_target.second) * x_est_target.lambda;

  auto gnss_position = main_observation.observed_values_.gnss_satellites_position.at(main_index);

  GnssObserveModel& main_observe_model = gnss_observed_models_.at(0);
  GnssObserveModel& target_observe_model = gnss_observed_models_.at(1);

  // SDCP
  z(row_offset) = carrier_phase_range_target - carrier_phase_range_main;
  h_x(row_offset) = target_observe_model.carrier_phase_range_model.at(target_index) - main_observe_model.carrier_phase_range_model.at(main_index);
  for (int j = 0; j < 3; ++j) {
    // position
    // ここの観測方程式が違う気がするな．
    const double los_element = (gnss_position[j] - x_est_main.position(j)) / main_observe_model.geometric_range.at(main_index); // 視線方向ベクトル要素
    // main
    H(row_offset, j) = los_element;
    // 定式化的には視線方向ベクトルは共通にすべき？あんまり関係ない気もする．
    // H(row_offset, NUM_SINGLE_STATE_ALL + j) = (x_est_target.position(j) - gnss_position[j]) / target_observe_model.geometric_range.at(target_index);
    // target
    H(row_offset, NUM_SINGLE_STATE_ALL + j) = -los_element;
  }
  // clock
  H(row_offset, 3) = -1.0;
  H(row_offset, NUM_SINGLE_STATE_ALL + 3) = 1.0;
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
  // この時のcommonがあるものから並んでいるとうれしい．<- 順番に関してはいったんかんがえんくていいか．
  for (const int& id : all_observed_gnss_ids)
  {
    // if main
    auto it_main = std::find(main_info_.now_observed_gnss_sat_id.begin(), main_info_.now_observed_gnss_sat_id.end(), id);
    if (it_main != main_info_.now_observed_gnss_sat_id.end() && *it_main == id)
    {
      UpdateObservationsGRAPHIC(0, x_est_main, id, z, h_x, H, Rv);
    }
    // if target
    auto it_target = std::find(target_info_.now_observed_gnss_sat_id.begin(), target_info_.now_observed_gnss_sat_id.end(), id);
    if (it_target != target_info_.now_observed_gnss_sat_id.end() && *it_target == id)
    {
      UpdateObservationsGRAPHIC(1, x_est_target, id, z, h_x, H, Rv);
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

// これは単にログ用のデータを更新しているだけ．ここでいいか．
void PBD_dgps::UpdateTrueAmbiguity(std::vector<std::vector<double>> N, const int gnss_sat_id, const double lambda)
{
  const int index_main = conv_index_from_gnss_sat_id(gnss_observations_.at(0).info_.now_observed_gnss_sat_id, gnss_sat_id);
  // ここは不定性だけを入れるように修正する．
  true_N_main(index_main) = N[0].at(gnss_sat_id);
  const int index_target = conv_index_from_gnss_sat_id(gnss_observations_.at(1).info_.now_observed_gnss_sat_id, gnss_sat_id);
  true_N_target(index_target) = N[1].at(gnss_sat_id);
  /*
  if (now_main_observing_ch.count(gnss_sat_id))
  {
    int ch = now_main_observing_ch[gnss_sat_id];
    true_N_main(ch) = N[0].at(gnss_sat_id)*lambda;
  }
  if (now_common_observing_ch.count(gnss_sat_id))
  {
    int ch = now_common_observing_ch[gnss_sat_id];
    true_N_target(ch) = N[1].at(gnss_sat_id)*lambda;
  }
  */
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

  const GnssObserveInfo& main_info_ = gnss_observations_.at(main_sat_id).info_;
  const GnssObserveInfo& target_info_ = gnss_observations_.at(target_sat_id).info_;
  // ここはiterで取得でいいのかも？
  for (int i = 0; i < main_info_.now_observed_gnss_sat_id.size(); ++i)
  {
    for (int j = 0; j < target_info_.now_observed_gnss_sat_id.size(); ++j)
    {
      int gnss_sat_id = main_info_.now_observed_gnss_sat_id.at(i);
      // どっかで複数衛星にも拡張
      if (gnss_sat_id == target_info_.now_observed_gnss_sat_id.at(j))
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
      //RemoveFromCh(gnss_sat_id, now_common_observing_ch, common_free_ch);
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

Eigen::VectorXd PBD_dgps::ConvReceivePosToCenterOfMass(Eigen::VectorXd x_state)
{
  Eigen::Vector3d pos_main = x_state.topRows(3);
  Quaternion q_i2b = main_dynamics_.GetQuaternion_i2b();
  libra::Vector<3> sat2ant_i = q_i2b.frame_conv_inv(antenna_pos_b_.at(0));
  for (uint8_t i = 0; i < 3; i++) pos_main(i) -= sat2ant_i[i]; // 補正する．

  Eigen::Vector3d pos_target = x_state.block(NUM_SINGLE_STATE_ALL, 0, 3, 1);
  q_i2b = target_dynamics_.GetQuaternion_i2b();
  sat2ant_i = q_i2b.frame_conv_inv(antenna_pos_b_.at(1));
  for (uint8_t i = 0; i < 3; i++) pos_target(i) -= sat2ant_i[i]; // 補正する．

  x_state.topRows(3) = pos_main;
  x_state.block(NUM_SINGLE_STATE_ALL, 0, 3, 1) = pos_target;

  return x_state;
}

Eigen::VectorXd PBD_dgps::ConvCenterOfMassToReceivePos(Eigen::VectorXd x_state)
{
  Eigen::Vector3d receive_pos_main = x_state.topRows(3);
  Quaternion q_i2b = main_dynamics_.GetQuaternion_i2b();
  libra::Vector<3> sat2ant_i = q_i2b.frame_conv_inv(antenna_pos_b_.at(0));
  for (uint8_t i = 0; i < 3; i++) receive_pos_main(i) += sat2ant_i[i]; // 補正する．

  Eigen::Vector3d receive_pos_target = x_state.block(NUM_SINGLE_STATE_ALL, 0, 3, 1);
  q_i2b = target_dynamics_.GetQuaternion_i2b();
  sat2ant_i = q_i2b.frame_conv_inv(antenna_pos_b_.at(1));
  for (uint8_t i = 0; i < 3; i++) receive_pos_target(i) += sat2ant_i[i]; // 補正する．

  x_state.topRows(3) = receive_pos_main;
  x_state.block(NUM_SINGLE_STATE_ALL, 0, 3, 1) = receive_pos_target;

  return x_state;
}


void PBD_dgps::UpdateBiasForm(const int sat_id, EstimatedVariables& x_est)// LEO衛星の数が増えたときは衛星ごとにこのクラスのインスタンスが生成される？ので一旦これで行く
{
  const GnssObserveInfo& observe_info_ = gnss_observations_.at(sat_id).info_;
  //観測する衛星同じだったら飛ばしていい
  if (CheckVectorEqual(observe_info_.pre_observed_gnss_sat_id, observe_info_.now_observed_gnss_sat_id))
  {
    return;
  }
  int n = observe_info_.now_observed_gnss_sat_id.size();
  int n_pre = observe_info_.pre_observed_gnss_sat_id.size();

  // index is the order in matrix
  int pre_index = 0;
  int now_index = 0;

  const std::vector<double> pre_estimated_bias = x_est.ambiguity.N; // 参照渡し？
  Eigen::MatrixXd pre_M = M;
  // reset
  M.block(sat_id*NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE, sat_id*NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE, NUM_GNSS_CH, NUM_GNSS_CH) = Eigen::MatrixXd::Zero(NUM_GNSS_CH, NUM_GNSS_CH);

  /*
  double geo_ure_sigma = sqrt(pow(L1_frequency / L2_frequency, 4.0) + 1.0) * pseudo_sigma / (pow(L1_frequency / L2_frequency, 2.0) - 1.0); //ionfree_pseudo_sigma
  */
  // i = gnss_sat_id
  for (int i = 0; i < num_of_gnss_satellites_; ++i)
  {
    if (observe_info_.pre_observed_status.at(i) == false && observe_info_.now_observed_status.at(i) == false) continue;
    // 見えなくなったとき
    else if (observe_info_.pre_observed_status.at(i) == true && observe_info_.now_observed_status.at(i) == false)
    {
      // 何もせず飛ばす．
      ++pre_index;
    }
    else if (observe_info_.pre_observed_status.at(i) == false && observe_info_.now_observed_status.at(i) == true)
    {
      // if (now_index != now_index_from_sat_it.at(i))
      // {
      //   cout << "now index something is wrong 1" << std::endl;
      //   abort();
      // }
      std::normal_distribution<> N_dist(0.0, sigma_N_ini);
      x_est.ambiguity.N.at(now_index) = N_dist(mt); // こいつもGauss Makov過程にした方がいいらしい?
      int offset = sat_id * NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE + now_index;

      // 行と列に関してもリセット
      M.block(0, offset, NUM_STATE_ALL, 1) = Eigen::MatrixXd::Zero(NUM_STATE_ALL, 1);
      M.block(offset, 0, 1, NUM_STATE_ALL) = Eigen::MatrixXd::Zero(1, NUM_STATE_ALL);
      M(offset, offset) = pow(sigma_N_ini, 2.0);
      ++now_index;
    }
    // 引き継ぐ
    else if (observe_info_.pre_observed_status.at(i) == true && observe_info_.now_observed_status.at(i) == true)
    {
      x_est.ambiguity.N.at(now_index) = pre_estimated_bias.at(pre_index);
      int offset_base = sat_id * NUM_SINGLE_STATE_ALL + NUM_SINGLE_STATE;
      M.block(0, offset_base + now_index, NUM_STATE_ALL, 1) = pre_M.block(0, offset_base + pre_index, NUM_STATE_ALL, 1);
      M.block(offset_base + now_index, 0, 1, NUM_STATE_ALL) = pre_M.block(offset_base + pre_index, 0, 1, NUM_STATE_ALL);
      ++pre_index; ++now_index;
    }
    if (now_index >= NUM_GNSS_CH || pre_index >= NUM_GNSS_CH) break; // ch以上の受信は出来ない
  }
  // now_index以上の部分を0に落とすということをやる．
  for (int i = now_index; i<NUM_GNSS_CH; ++i)
  x_est.ambiguity.N.at(i) = 0;
}

// これは観測情報を行列に入れている部分なので推定のところでするべき．
void PBD_dgps::SetBiasToObservation(const int sat_id, EstimatedVariables& x_est, PBD_GnssObservation& gnss_observation)
{
    // 共通衛星見つける
    FindCommonObservedGnss(std::make_pair(0, 1)); // ここではない気がする．

    Eigen::MatrixXd pre_M = M;
    int n = gnss_observation.info_.now_observed_gnss_sat_id.size();
    int n_pre = gnss_observation.info_.pre_observed_gnss_sat_id.size();

    UpdateBiasForm(sat_id, x_est);
    gnss_observation.UpdateInfoAfterObserved();

    return;
}

//使ってない，修正してこの形に持っていく
Eigen::VectorXd PBD_dgps::CalculateSingleDifference(const Eigen::VectorXd& main_observation, const Eigen::VectorXd& target_observation) const
{
  return main_observation - target_observation;
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

void PBD_dgps::DynamicNoiseScaling(Eigen::MatrixXd Q_dash, Eigen::MatrixXd Phi, Eigen::MatrixXd H)
{
  // 観測更新後のPに対して行う．
  Eigen::MatrixXd P_dash = Phi * M * Phi.transpose() + Q_dash;
  Eigen::MatrixXd P      = Phi * M * Phi.transpose() + Q;
#define TRACE_SCALE_

#ifdef TRACE_SCALE_
  double beta_dash = (H * P_dash * H.transpose()).trace() / (H * P * H.transpose()).trace();
  Q = sqrt(beta_dash) * Q;
#else
  // traceとらない方法.位置精度に依存する．
  Eigen::VectorXd diag_dash = P_dash.diagonal();
  Eigen::VectorXd diag = P.diagonal();
  for (uint8_t i = 0; i < diag.size(); i++)
  {
    if (diag(i) == 0) continue;
    Q(i, i) = sqrt(diag_dash(i) / diag(i)) * Q(i, i);
  }
#endif
}
