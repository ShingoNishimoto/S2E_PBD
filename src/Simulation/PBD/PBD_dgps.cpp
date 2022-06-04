#include <iomanip>
#include <set>
#include "PBD_dgps.h"
#include "PBD_const.h"

// clock model
#define CLOCK_IS_WHITE_NOISE (1)
// #define CLOCK_IS_RANDOM_WALK (1)

// inverse calculation method
#define CHOLESKY (0)
#define QR (1)

static const int conv_index_from_gnss_sat_id(std::vector<int> observed_gnss_sat_id, const int gnss_sat_id);

// outputを変えるときは"result.csv"を変更する．せめてパスは変えたい．
PBD_dgps::PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Orbit& main_orbit, const Orbit& target_orbit, PBD_GnssObservation& main_observation, PBD_GnssObservation& target_observation) :mt(42), step_time(sim_time_.GetStepSec()), ofs("result_new.csv"), num_of_gnss_satellites_(gnss_satellites_.GetNumOfSatellites()), main_orbit_(main_orbit), target_orbit_(target_orbit), receiver_clock_bias_main_(main_observation.receiver_clock_bias_), receiver_clock_bias_target_(target_observation.receiver_clock_bias_), gnss_observations_({main_observation, target_observation})
{
  //初期化
  x_est_main.position = Eigen::VectorXd::Zero(3);
  libra::Vector<3> position_main = main_orbit_.GetSatPosition_i();
  for (int i = 0; i < 3; ++i) x_est_main.position(i) = position_main[i];
  x_est_main.clock = Eigen::VectorXd::Zero(1);
  libra::Vector<3> velocity_main = main_orbit_.GetSatVelocity_i();
  for (int i = 0; i < 3; ++i) x_est_main.velocity(i) = velocity_main[i];
  x_est_main.acceleration = Eigen::VectorXd::Zero(3);
  x_est_main.bias = Eigen::VectorXd::Zero(num_of_gnss_channel);

  x_est_target.position = Eigen::VectorXd::Zero(3);
  libra::Vector<3> position_target = target_orbit_.GetSatPosition_i();
  for (int i = 0; i < 3; ++i) x_est_target.position(i) = position_target[i];
  x_est_target.clock = Eigen::VectorXd::Zero(1);
  libra::Vector<3> velocity_target = target_orbit_.GetSatVelocity_i();
  for (int i = 0; i < 3; ++i) x_est_target.velocity(i) = velocity_target[i];
  x_est_target.acceleration = Eigen::VectorXd::Zero(3);
  x_est_target.bias = Eigen::VectorXd::Zero(num_of_gnss_channel);

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
  
  true_bias_main = Eigen::VectorXd::Zero(num_of_gnss_channel);
  true_bias_target = Eigen::VectorXd::Zero(num_of_gnss_channel);

  // 初期分散
  std::normal_distribution<> position_dist(0.0,sigma_r_ini);
  std::normal_distribution<> receiver_clock_dist(0.0, sigma_cdt_ini);
  std::normal_distribution<> velocity_dist(0.0, sigma_v_ini);
  std::normal_distribution<> acc_r_dist(0.0, sigma_acc_r_ini);
  std::normal_distribution<> acc_t_dist(0.0, sigma_acc_t_ini);
  std::normal_distribution<> acc_n_dist(0.0, sigma_acc_n_ini);
  std::normal_distribution<> N_dist(0.0, sigma_N_ini);

  // Vって名前は変えた方がいいかも
  Eigen::VectorXd V = Eigen::VectorXd::Constant(state_dimension, 0);

  // A-priori 
  for(int i = 0; i < 3; ++i) V(i) = pow(sigma_r_ini, 2.0); // position
  V(3) = pow(sigma_cdt_ini, 2.0); // clock
  for(int i = 0; i < 3; ++i) V(4 + i) = pow(sigma_v_ini, 2.0); // velocity
  V(7) = pow(sigma_acc_r_ini, 2.0);
  V(8) = pow(sigma_acc_t_ini, 2.0);
  V(9) = pow(sigma_acc_n_ini, 2.0);
  // 初期は全部に入れている．
  for (int i = 0; i < num_of_gnss_channel; ++i) V(num_of_single_status + i) = pow(sigma_N_ini, 2.0); // N
  // 以下がtarget
  for (int i = 0; i < 3; ++i) V(single_dimension + i) = pow(sigma_r_ini, 2.0);
  V(single_dimension + 3) = pow(sigma_cdt_ini, 2.0);
  for (int i = 0; i < 3; ++i) V(single_dimension + 4 + i) = pow(sigma_v_ini, 2.0);
  V(single_dimension + 7) = pow(sigma_acc_r_ini, 2.0);
  V(single_dimension + 8) = pow(sigma_acc_t_ini, 2.0);
  V(single_dimension + 9) = pow(sigma_acc_n_ini, 2.0);
  for (int i = 0; i < num_of_gnss_channel; ++i) V(single_dimension + num_of_single_status + i) = pow(sigma_N_ini, 2.0); // N

  M = V.asDiagonal(); // 誤差分散共分散行列M

  // 初期位置はガウシアンからサンプリング．mtは乱数のシード
  for(int i = 0; i < 3; ++i) x_est_main.position(i) += position_dist(mt);
  x_est_main.clock(0) += receiver_clock_dist(mt);
  for (int i = 0; i < 3; ++i) x_est_main.velocity(i) += velocity_dist(mt);
  // for(int i = 0; i < 3; ++i) x_est_main.acceleration(i) += acc_r_dist(mt);
  x_est_main.acceleration(0) += acc_r_dist(mt);
  x_est_main.acceleration(1) += acc_t_dist(mt);
  x_est_main.acceleration(2) += acc_n_dist(mt);
  for(int i = 0; i < num_of_gnss_channel; ++i) x_est_main.bias(i) += N_dist(mt);
  for(int i = 0; i < 3; ++i) x_est_target.position(i) += position_dist(mt);
  x_est_target.clock(0) += receiver_clock_dist(mt);
  for (int i = 0; i < 3; ++i) x_est_target.velocity(i) += velocity_dist(mt);
  // for(int i = 0; i < 3; ++i) x_est_target.acceleration(i) += acc_r_dist(mt);
  x_est_target.acceleration(0) += acc_r_dist(mt);
  x_est_target.acceleration(1) += acc_t_dist(mt);
  x_est_target.acceleration(2) += acc_n_dist(mt);
  for(int i = 0; i < num_of_gnss_channel; ++i) x_est_target.bias(i) += N_dist(mt);

  common_observed_status.assign(num_of_gnss_satellites_, false);
  // for (int i = 0; i < num_of_gnss_channel; ++i) main_free_ch.push_back(i);
  // for (int i = 0; i < num_of_gnss_channel; ++i) common_free_ch.push_back(i);

  std::ofstream ofs_ini_txt("readme_new.txt");
  ofs_ini_txt << "initial position dist: " << sigma_r_ini << std::endl;
  ofs_ini_txt << "initial velocity dist: " << sigma_v_ini << std::endl;
  ofs_ini_txt << "initial acceleration dist: " << sigma_acc_r_ini << std::endl;
  ofs_ini_txt << "initial clock dist: " << sigma_cdt_ini << std::endl;
  ofs_ini_txt << "initial ambiguity dist: " << sigma_N_ini << std::endl;
  ofs_ini_txt << "pseudo dist: " << pseudo_sigma << std::endl;
  ofs_ini_txt << "carrier dist: " << carrier_sigma << std::endl;
  ofs_ini_txt << "clock dist: " << clock_sigma << std::endl;
  ofs_ini_txt << "process noise of position: " << sigma_r_process << std::endl;
  ofs_ini_txt << "process noise of velocity: " << sigma_v_process << std::endl;
  ofs_ini_txt << "process noise of radial acceleration: " << sigma_acc_r_process << std::endl;
  ofs_ini_txt << "process noise of tangential acceleration: " << sigma_acc_t_process << std::endl;
  ofs_ini_txt << "process noise of north acceleration: " << sigma_acc_n_process << std::endl;
  ofs_ini_txt << "process noise of clock: " << sigma_cdt_process << std::endl;
  ofs_ini_txt << "process noise of ambiguity: " << sigma_N_process << std::endl;
  ofs_ini_txt << "time const. acceleration: " << tau_a << std::endl;
  ofs_ini_txt << "time const. clock: " << tau_cdt << std::endl;
  ofs_ini_txt << "mask angle: " << gnss_observations_.at(0).mask_angle << std::endl; // FIXME
  ofs_ini_txt << "num of status: " << num_of_status << std::endl;
  ofs_ini_txt << "observe step time: " << observe_step_time << std::endl;
  ofs_ini_txt << "log step time: " << log_step_time << std::endl;
}

PBD_dgps::~PBD_dgps(){}

void PBD_dgps::Update(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, PBD_GnssObservation& main_observation_, PBD_GnssObservation& target_observation_)// , const Orbit& main_orbit, const Orbit& target_orbit)
{
  // 参照渡ししたものを代入するのは無理．
  gnss_observations_.clear();
  gnss_observations_.push_back(main_observation_);
  gnss_observations_.push_back(target_observation_);

  double elapsed_time = sim_time_.GetElapsedSec();
  double tmp = floor(elapsed_time/observe_step_time + 1e-4); //1e-4は数値誤差
  double tmp_log = floor(elapsed_time/log_step_time + 1e-4);

  //まず更新
  OrbitPropagation();

  //観測時間にピッタリ
  if(abs(elapsed_time - tmp*observe_step_time) < step_time/2.0){

    // ここの引数をどう渡すかとかが関係している？
    // PBD_GnssObservation& main_observation_ = gnss_observations_.at(0);
    SetBiasToObservation(0, x_est_main, gnss_observations_.at(0));
    // PBD_GnssObservation& target_observation_ = gnss_observations_.at(1);
    SetBiasToObservation(1, x_est_target, gnss_observations_.at(1));
    // ここでinfoの更新がしたい．
    KalmanFilter();
  }

  //log output
  if(abs(elapsed_time - tmp_log*log_step_time) < step_time/2.0){
    libra::Vector<3> sat_position = main_orbit_.GetSatPosition_i();
    libra::Vector<3> sat_velocity = main_orbit_.GetSatVelocity_i();
    for(int i = 0;i < 3;++i) ofs << std::fixed << std::setprecision(30) << sat_position[i] << ","; // r_m_true
    ofs << std::fixed << std::setprecision(30) << receiver_clock_bias_main_ << ","; // t_m_true
    for(int i = 0;i < 3;++i) ofs << std::fixed << std::setprecision(30) << sat_velocity[i] << ","; // v_m_ture
    libra::Vector<3> sat_position_target = target_orbit_.GetSatPosition_i();
    libra::Vector<3> sat_velocity_target = target_orbit_.GetSatVelocity_i();
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(30) << sat_position_target[i] << ","; // r_t_ture
    ofs << std::fixed << std::setprecision(30) << receiver_clock_bias_target_ << ","; // t_t_true
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(30) << sat_velocity_target[i] << ","; // v_t_ture
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(30) << x_est_main.position(i) << ","; // r_m_est
    ofs << std::fixed << std::setprecision(30) << x_est_main.clock(0) << ","; // t_m_est
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(30) << x_est_main.velocity(i) << ","; // v_m_est
    for(int i = 0;i < 3;++i) ofs << std::fixed << std::setprecision(30) << x_est_main.acceleration(i) << ","; // a_m_est
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(30) << x_est_target.position(i) << ","; // r_t_est
    ofs << std::fixed << std::setprecision(30) << x_est_target.clock(0) << ","; //t_t_est
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(30) << x_est_target.velocity(i) << ","; // v_t_est
    for (int i = 0; i < 3; ++i) ofs << std::fixed << std::setprecision(30) << x_est_target.acceleration(i) << ","; // a_t_est
    for (int i = 0; i < num_of_gnss_channel; ++i) ofs << std::fixed << std::setprecision(30) << true_bias_main(i) << ","; // N_true
    for (int i = 0; i < num_of_gnss_channel; ++i) ofs << std::fixed << std::setprecision(30) << x_est_main.bias(i) << ","; // N_est
    for (int i = 0; i < num_of_gnss_channel; ++i) ofs << std::fixed << std::setprecision(30) << true_bias_target(i) << ","; // N_true
    for (int i = 0; i < num_of_gnss_channel; ++i) ofs << std::fixed << std::setprecision(30) << x_est_target.bias(i) << ","; // N_est
    for (int i = 0; i < state_dimension; ++i) ofs << std::fixed << std::setprecision(30) << M(i, i) << ",";
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
    for (int i = 0; i < num_of_gnss_channel; ++i)
    {
      if (i >= visible_gnss_num_main) ofs << -1 << ",";
      else ofs << gnss_observations_.at(0).info_.now_observed_gnss_sat_id.at(i) << ",";
    }
    // target observe gnss sat id
    for (int i = 0; i < num_of_gnss_channel; ++i)
    {
      if (i >= visible_gnss_num_target) ofs << -1 << ",";
      else ofs << gnss_observations_.at(1).info_.now_observed_gnss_sat_id.at(i) << ",";
    }
    ofs << std::endl;
  }

  return;
}

void PBD_dgps::OrbitPropagation()
{
  //RK4
  Eigen::Vector3d position_main = x_est_main.position;
  Eigen::Vector3d velocity_main = x_est_main.velocity;
  Eigen::Vector3d acceleration_main = x_est_main.acceleration;
  Eigen::Vector3d position_target = x_est_target.position;
  Eigen::Vector3d velocity_target = x_est_target.velocity;
  Eigen::Vector3d acceleration_target = x_est_target.acceleration;

  // ここはcoreの機能を使うように修正
  std::vector<Eigen::Vector3d> pos_vel_main = RK4(position_main, velocity_main, acceleration_main);
  position_main = pos_vel_main[0];
  velocity_main = pos_vel_main[1];

  std::vector<Eigen::Vector3d> pos_vel_target = RK4(position_target, velocity_target, acceleration_target);
  position_target = pos_vel_target[0];
  velocity_target = pos_vel_target[1];

  x_est_main.position = position_main;
  x_est_main.velocity = velocity_main;
  x_est_target.position = position_target;
  x_est_target.velocity = velocity_target;

  // acceleration ここでは q = sigma_acc_process
  Eigen::MatrixXd Phi_a = CalculatePhi_a(step_time);
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
  // cdt
// #ifdef CLOCK_IS_RANDOM_WALK
//   std::normal_distribution<> cdt_process_noise(0.0, sigma_cdt_process*sqrt(step_time/tau_cdt));
//   std::normal_distribution<> cdt_process_noise(0.0, sigma_cdt_process);
//   x_est_main.clock(0) = cdt_process_noise(mt);
//   x_est_target.clock(0) = cdt_process_noise(mt);
// #endif // CLOCK_IS_RANDOM_WALK
  M = UpdateM();
}

// Orbitの更新を使えるように修正．位置と速度はRK4で伝搬しているのが，共分散はオイラー法になっている．EKFなので仕方がない．
std::vector<Eigen::Vector3d> PBD_dgps::RK4(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration)
{
  Eigen::Vector3d k0 = PositionDifferential(velocity);
  Eigen::Vector3d l0 = VelocityDifferential(position, velocity, acceleration);

  Eigen::Vector3d tmp_position = position + k0 * step_time / 2.0;
  Eigen::Vector3d tmp_velocity = velocity + l0 * step_time / 2.0;
  Eigen::Vector3d k1 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l1 = VelocityDifferential(tmp_position, tmp_velocity, acceleration);

  tmp_position = position + k1 * step_time / 2.0;
  tmp_velocity = velocity + l1 * step_time / 2.0;
  Eigen::Vector3d k2 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l2 = VelocityDifferential(tmp_position, tmp_velocity, acceleration);

  tmp_position = position + k2 * step_time;
  tmp_velocity = velocity + l2 * step_time;
  Eigen::Vector3d k3 = PositionDifferential(tmp_velocity);
  Eigen::Vector3d l3 = VelocityDifferential(tmp_position, tmp_velocity, acceleration);

  std::vector<Eigen::Vector3d> position_velocity = 
  { position + step_time * (k0 + 2.0 * k1 + 2.0 * k2 + k3) / 6.0, 
    velocity + step_time * (l0 + 2.0 * l1 + 2.0 * l2 + l3) / 6.0 };
  return position_velocity;
}

Eigen::Vector3d PBD_dgps::PositionDifferential(const Eigen::Vector3d& velocity) const
{
  return velocity;
}

// 単位が違うのでそこの変換が必要　
Eigen::Vector3d PBD_dgps::VelocityDifferential(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration) const
{
  double r = position.norm();
  double v = velocity.norm();

  double z = position(2);
  double conv_nm2m = 1e-6; // accの単位補正

  double ac_norm = - mu_const/position.squaredNorm(); //2体の重力項
  double tmp_J2_coefficient = 3.0/2.0*mu_const*J2_const*pow(Earth_Radius, 2.0)/pow(r, 4.0); //J2項の係数

  Eigen::Vector3d all_acceleration = position/r;

  all_acceleration(0) *= ac_norm - tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0));
  all_acceleration(1) *= ac_norm - tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0));
  all_acceleration(2) *= ac_norm - tmp_J2_coefficient*(3.0 - 5.0*pow(z/r, 2.0));

  all_acceleration -= Cd*v*velocity; //-Cd*V^2*(Vi/V) 大気抵抗

  // ここも変換が必要．
  
  Eigen::MatrixXd acc(3, 1);
  acc.block(0, 0, 3, 1) = acceleration;
  Eigen::Vector3d acc_eci = TransRTN2ECI(position, velocity)*acc;
  all_acceleration +=  acc_eci*conv_nm2m; //残りの摂動要素
  // all_acceleration +=  acceleration * conv_nm2m; //残りの摂動要素

  return all_acceleration; // m/s2
}

Eigen::MatrixXd PBD_dgps::UpdateM()
{
  // int n_main = estimated_bias.size(); //mainの見るGNSS数
  // int n_common = estimated_bias.size(); //共通のGNSS数

  // Dinamics and Kinematics model
  Eigen::MatrixXd A_jacobi = CalculateA(x_est_main, x_est_target); //Jacobi行列
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(state_dimension, state_dimension);
  A.topLeftCorner(num_of_single_status, num_of_single_status) = A_jacobi.topLeftCorner(num_of_single_status, num_of_single_status);
  A.block(single_dimension, single_dimension, num_of_single_status, num_of_single_status) = A_jacobi.bottomRightCorner(num_of_single_status, num_of_single_status); 
  // STM
  Eigen::MatrixXd Phi = Eigen::MatrixXd::Identity(state_dimension, state_dimension); // Nのところは単位行列
  Phi += step_time * A;
  // clock
#ifdef CLOCK_IS_RANDOM_WALK
  Phi(3, 3) = 1.0; // random walkなのでΦは1
  Phi(single_dimension + 3, single_dimension + 3) = 1.0;
#endif // CLOCK_IS_RANDOM_WALK
#ifdef CLOCK_IS_WHITE_NOISE
  Phi(3, 3) = 0; // white noise model
  Phi(single_dimension + 3, single_dimension + 3) = 0;
#endif // CLOCK_IS_WHITE_NOISE
  // a
  Phi.block(7, 7, 3, 3) = CalculatePhi_a(step_time);
  Phi.block(single_dimension + 7, single_dimension + 7, 3, 3) = CalculatePhi_a(step_time);

  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(state_dimension, 8); // aとclock
  B(3, 0) = 1.0;
  B(single_dimension +  3, 4) = 1.0;
  B.block(7, 1, 3, 3) = Eigen::Matrix3d::Identity();
  B.block(single_dimension + 7, 5, 3, 3) = Eigen::Matrix3d::Identity();

  // ここもaについて修正．
  Eigen::MatrixXd Q = CalculateQ();
  // 観測するたびにNの部分を初期化？

  // Gamma = BQB^t
  Eigen::MatrixXd Gamma = B * Q * B.transpose();
  // add process noise for r and v
  Gamma.block(0, 0, 3, 3) = pow(sigma_r_process, 2.0) * Eigen::Matrix3d::Identity()* pow(step_time, 2.0);
  Gamma.block(4, 4, 3, 3) = pow(sigma_v_process, 2.0) * Eigen::Matrix3d::Identity()* pow(step_time, 2.0);
  Gamma.block(single_dimension, single_dimension, 3, 3) = pow(sigma_r_process, 2.0) * Eigen::Matrix3d::Identity()*pow(step_time, 2.0);
  Gamma.block(single_dimension + 4, single_dimension+ 4, 3, 3) = pow(sigma_v_process, 2.0) * Eigen::Matrix3d::Identity()* pow(step_time, 2.0);
  // N process <- ここやる！！！！

  Eigen::MatrixXd res = Phi * M * Phi.transpose() + Gamma;
  // Nのない部分を0に落とす
  int n_main = gnss_observations_.at(0).info_.now_observed_gnss_sat_id.size();
  res.block(num_of_single_status + n_main, num_of_single_status + n_main, num_of_gnss_channel - n_main, num_of_gnss_channel - n_main) = Eigen::MatrixXd::Zero(num_of_gnss_channel - n_main, num_of_gnss_channel - n_main);
  int n_target = gnss_observations_.at(1).info_.now_observed_gnss_sat_id.size();
  res.block(single_dimension + num_of_single_status + n_target, single_dimension + num_of_single_status + n_target, num_of_gnss_channel - n_target, num_of_gnss_channel - n_target) = Eigen::MatrixXd::Zero(num_of_gnss_channel - n_target, num_of_gnss_channel - n_target);
  
  return res;
}

// この段階必要ない気がする．
Eigen::MatrixXd PBD_dgps::CalculateA(const EstimatedVariables& x_est_main, const EstimatedVariables& x_est_target) const
{

  Eigen::MatrixXd A_main = CalculateJacobian(x_est_main.position, x_est_main.velocity, x_est_main.acceleration);
  Eigen::MatrixXd A_target = CalculateJacobian(x_est_target.position, x_est_target.velocity, x_est_target.acceleration);
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_of_status, num_of_status);
  A.topLeftCorner(num_of_single_status, num_of_single_status) = A_main;
  A.bottomRightCorner(num_of_single_status, num_of_single_status) = A_target;

  // Kinematics model 多分これやと誤差大きすぎて無理な気がする．同じようにDynamics使う必要がある．
  // (dr, dv)
  // A(num_of_single_status, num_of_single_status + 4)     = 1.0;
  // A(num_of_single_status + 1, num_of_single_status + 5) = 1.0;
  // A(num_of_single_status + 2, num_of_single_status + 6) = 1.0;
  // // (dv, da)
  // A(num_of_single_status + 4, num_of_single_status + 7) = 1.0 * 1e-6;
  // A(num_of_single_status + 5, num_of_single_status + 8) = 1.0 * 1e-6;
  // A(num_of_single_status + 6, num_of_single_status + 9) = 1.0 * 1e-6;
  // // (dr, da)
  // A(num_of_single_status, num_of_single_status + 7) = 1.0 * 1e-6* step_time/2;
  // A(num_of_single_status + 1, num_of_single_status + 8) = 1.0 * 1e-6* step_time/2;
  // A(num_of_single_status + 2, num_of_single_status + 9) = 1.0 * 1e-6* step_time/2;

  return A;
}

Eigen::MatrixXd PBD_dgps::CalculateJacobian(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration) const
{ 
  double r = position.norm(); // [m]
  double v = velocity.norm(); // [m/s]
  double a = acceleration.norm(); // [nm/s]

  double x = position(0); double y = position(1); double z = position(2);
  double vx = velocity(0); double vy = velocity(1); double vz = velocity(2);

  double J2_coefficient = 3.0 / 2.0 * mu_const * J2_const * Earth_Radius * Earth_Radius;

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_of_single_status, num_of_single_status);
  // (r, v)
  A(0, 4) = 1.0; A(1, 5) = 1.0; A(2, 6) = 1.0;
  // (v, r)
  A(4, 0) = 3.0 * mu_const * x * x / pow(r, 5.0) - mu_const / pow(r, 3.0) - J2_coefficient * (1.0 / pow(r, 5.0) - 5.0 * (x * x + z * z) / pow(r, 7.0) + 35.0 * x * x * z * z / pow(r, 9.0));
  A(4, 1) = 3.0 * mu_const * x * y / pow(r, 5.0) - J2_coefficient * (-5.0 * x * y / pow(r, 7.0) + 35.0 * x * y * z * z / pow(r, 9.0));
  A(4, 2) = 3.0 * mu_const * x * z / pow(r, 5.0) - J2_coefficient * (-15.0 * x * z / pow(r, 7.0) + 35.0 * x * z * z * z / pow(r, 9.0));
  A(5, 0) = 3.0 * mu_const * x * y / pow(r, 5.0) - J2_coefficient * (-5.0 * x * y / pow(r, 7.0) + 35.0 * x * y * z * z / pow(r, 9.0));
  A(5, 1) = 3.0 * mu_const * y * y / pow(r, 5.0) - mu_const / pow(r, 3.0) - J2_coefficient * (1.0 / pow(r, 5.0) - 5.0 * (y * y + z * z) / pow(r, 7.0) + 35.0 * y * y * z * z / pow(r, 9.0));
  A(5, 2) = 3.0 * mu_const * y * z / pow(r, 5.0) - J2_coefficient * (-15.0 * y * z / pow(r, 7.0) + 35.0 * y * z * z * z / pow(r, 9.0));
  A(6, 0) = 3.0 * mu_const * x * z / pow(r, 5.0) - J2_coefficient * (-15.0 * x * z / pow(r, 7.0) + 35.0 * x * z * z * z / pow(r, 9.0));
  A(6, 1) = 3.0 * mu_const * y * z / pow(r, 5.0) - J2_coefficient * (-15.0 * y * z / pow(r, 7.0) + 35.0 * y * z * z * z / pow(r, 9.0));
  A(6, 2) = 3.0 * mu_const * z * z / pow(r, 5.0) - mu_const / pow(r, 3.0) - J2_coefficient * (3.0 / pow(r, 5.0) - 30.0 * z * z / pow(r, 7.0) + 35.0 * pow(z, 4.0) / pow(r, 9.0));

  // これ速度に入る？普通は加速度として入ってくるはず
  // A(4, 4) = -Cd * (vx * vx / v + v);    A(4, 5) = -Cd * vx * vy / v;    A(4, 6) = -Cd * vx * vz / v;
  // A (5, 4) = -Cd * vx * vy / v;    A(5, 5) = -Cd * (vy * vy / v + v);    A(5, 6) = -Cd * vy * vz / v;
  // A(6, 4) = -Cd * vx * vz / v;    A(6, 5) = -Cd * vy * vz / v;    A(6, 6) = -Cd * (vz * vz / v + v);

  // (v, a)
  Eigen::Matrix3d rtn2eci = TransRTN2ECI(position, velocity);
  A.block(4, 7, 3, 3) = 1e-6*rtn2eci;
  // A(4, 7) = 1.0*1e-6;	A(5, 8) = 1.0 * 1e-6;	A(6, 9) = 1.0 * 1e-6;

  // これはもしCdも推定しているならいる．
  //A(4,10) = -v*vx;    A(5,10) = -v*vy;    A(6,10) = -v*vz;
  return A;
};

Eigen::Matrix3d PBD_dgps::TransRTN2ECI(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const
{
  Eigen::MatrixXd rtn_r(3, 1);
  for (int i = 0; i < 3;++i) rtn_r(i, 0) = position(i) / position.norm();
  Eigen::MatrixXd rtn_t(3, 1);
  for (int i = 0; i < 3; ++i) rtn_t(i, 0) = velocity(i) / velocity.norm();
  // Eigen::Vector3d rtn_n = rtn_r.cross(rtn_t); // cross使えんのクソ
  Eigen::MatrixXd rtn_n(3, 1);
  rtn_n(0, 0) = rtn_r(1, 0) * rtn_t(2, 0) - rtn_r(2, 0) * rtn_t(1, 0);
  rtn_n(1, 0) = rtn_r(2, 0) * rtn_t(0, 0) - rtn_r(0, 0) * rtn_t(2, 0);
  rtn_n(2, 0) = rtn_r(0, 0) * rtn_t(1, 0) - rtn_r(1, 0) * rtn_t(0, 0);
  rtn_n /= rtn_n.norm();
  Eigen::MatrixXd RTN2ECI(3,3);
  RTN2ECI.block(0, 0, 3, 1) = rtn_r;
  RTN2ECI.block(0, 1, 3, 1) = rtn_t;
  RTN2ECI.block(0, 2, 3, 1) = rtn_n;
  return RTN2ECI;
};

// Process noiseのvarianceを計算．Bを使う形に修正．
Eigen::MatrixXd PBD_dgps::CalculateQ()
{
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(8, 8);
  double phi = exp(-step_time / tau_a);
  double q_acc_r = pow(sigma_acc_r_process, 2.0) * (1 - pow(phi, 2.0));
  double q_acc_t = pow(sigma_acc_t_process, 2.0) * (1 - pow(phi, 2.0));
  double q_acc_n = pow(sigma_acc_n_process, 2.0) * (1 - pow(phi, 2.0));
  double q_cdt;
#ifdef CLOCK_IS_WHITE_NOISE
  q_cdt = pow(sigma_cdt_process, 2.0);
#endif // CLOCK_IS_WHITE_NOISE
#ifdef CLOCK_IS_RANDOM_WALK
  q_cdt = pow(sigma_cdt_process, 2.0) * (step_time / tau_cdt);
#endif // CLOCK_IS_RANDOM_WALK

  // for(int i = 0;i < 3;++i) Q(i, i) = pow(sigma_r_process, 2.0);
  Q(0, 0) = q_cdt;
  // for(int i = 4;i < 7;++i) Q(i, i) = pow(sigma_v_process, 2.0);
  Q(1, 1) = q_acc_r;
  Q(2, 2) = q_acc_t;
  Q(3, 3) = q_acc_n;
  // for (int i = num_of_single_status; i < n_main + num_of_single_status; ++i) Q(i, i) = pow(sigma_N_process, 2.0);
  // for (int i = num_of_single_status; i < num_of_single_status + 3; ++i) Q(i, i) = pow(sigma_r_process, 2.0);
  Q(4, 4) = q_cdt;
  // for (int i = num_of_single_status + 4; i < num_of_single_status + 7; ++i) Q(i, i) = pow(sigma_v_process, 2.0);
  Q(5, 5) = q_acc_r;
  Q(6, 6) = q_acc_t;
  Q(7, 7) = q_acc_n;
  // for (int i = num_of_single_status + single_offset; i < num_of_status + n_main + n_common; ++i) Q(i, i) = pow(sigma_N_process, 2.0);

  return Q;
}

Eigen::MatrixXd PBD_dgps::CalculatePhi_a(const double dt)
{
  // ここが各軸に応じて変わる部分なのか．
  double phi = exp(-dt / tau_a); // constやないか
  Eigen::MatrixXd Phi_a = Eigen::MatrixXd::Identity(3, 3);
  Phi_a *= phi;
  return Phi_a;
};

void PBD_dgps::KalmanFilter()
{
  // gnss_observed_hoge の観測される時間はどうなっているのか？
  // int n_common = common_observed_gnss_sat_id.size();
  
  // GRAPHIC*2 + SDCPにする．
  Eigen::VectorXd z = Eigen::VectorXd::Zero(observation_dimension); //観測ベクトル
  Eigen::VectorXd h_x = Eigen::VectorXd::Zero(observation_dimension); // 観測モデル行列

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(observation_dimension, state_dimension); //観測行列（dhx/dx）
  // Oneにしている部分がどう効いてくるのかの確認は必要．
  Eigen::VectorXd R_V = Eigen::MatrixXd::Zero(observation_dimension, 1); // 観測誤差共分散．
  const GnssObserveInfo& main_info_ = gnss_observations_.at(0).info_;
  const GnssObserveInfo& target_info_ = gnss_observations_.at(1).info_;
  std::vector<int> all_observed_gnss_ids = main_info_.now_observed_gnss_sat_id;
  all_observed_gnss_ids.insert(all_observed_gnss_ids.end(), target_info_.now_observed_gnss_sat_id.begin(), target_info_.now_observed_gnss_sat_id.end()); // concate
  sort(all_observed_gnss_ids.begin(), all_observed_gnss_ids.end());
  all_observed_gnss_ids.erase(unique(all_observed_gnss_ids.begin(), all_observed_gnss_ids.end()), all_observed_gnss_ids.end()); // unique
  // この時のcommonがあるものから並んでいるとうれしい．<- 順番に関してはいったんかんがえんくていいか．
  for (const int& id :all_observed_gnss_ids)
  {
    // if main
    auto it_main = std::find(main_info_.now_observed_gnss_sat_id.begin(), main_info_.now_observed_gnss_sat_id.end(), id);
    if (it_main != main_info_.now_observed_gnss_sat_id.end() && *it_main == id)
    {
      UpdateObservationsGRAPHIC(0, x_est_main, id, z, h_x, H, R_V);
    }
    // if target
    auto it_target = std::find(target_info_.now_observed_gnss_sat_id.begin(), target_info_.now_observed_gnss_sat_id.end(), id);
    if (it_target != target_info_.now_observed_gnss_sat_id.end() && *it_target == id)
    {
      UpdateObservationsGRAPHIC(1, x_est_target, id, z, h_x, H, R_V);
    }
    // if common
    auto it_common = std::find(common_observed_gnss_sat_id.begin(), common_observed_gnss_sat_id.end(), id);
    if (it_common != common_observed_gnss_sat_id.end() && *it_common == id)
    {
      UpdateObservationsSDCP(id, z, h_x, H, R_V);
    }
  }
    
  Eigen::MatrixXd R = R_V.asDiagonal();

  Eigen::MatrixXd hmh = H * M * H.transpose();
  /*
  if (abs(hmh.determinant()) < 10e-10)
  {
    cout << "HMHt matirx is singular" << std::endl;
    abort();
  }
  */
  Eigen::MatrixXd tmp = R + hmh; // (observation_num, observation_num)
  // 観測量がないところは0にし直すみたいな加工が必要かもしれない．数字の桁数で打ち切りみたいな形にできればいい．
  //カルマンゲイン
  Eigen::MatrixXd K = CalculateK(H, tmp);

  Eigen::VectorXd x_predict = Eigen::VectorXd(state_dimension);
  x_predict.topRows(3) = x_est_main.position;
  x_predict.block(3, 0, 1, 1) = x_est_main.clock;
  x_predict.block(4, 0, 3, 1)  = x_est_main.velocity;
  x_predict.block(7, 0, 3, 1) = x_est_main.acceleration;
  x_predict.block(10, 0, num_of_gnss_channel, 1) = x_est_main.bias;
  //x_predict(10) = Cd;
  x_predict.block(single_dimension, 0, 3, 1) = x_est_target.position;
  x_predict.block(single_dimension + 3, 0, 1, 1) = x_est_target.clock;
  x_predict.block(single_dimension + 4, 0, 3, 1) = x_est_target.velocity;
  x_predict.block(single_dimension + 7, 0, 3, 1) = x_est_target.acceleration;
  x_predict.bottomRows(num_of_gnss_channel) = x_est_target.bias;

  Eigen::VectorXd x_update = x_predict + K*(z - h_x); // 完全な0ではないからリセットが必要?

  //更新
  x_est_main.position = x_update.topRows(3);
  x_est_main.clock = x_update.block(3, 0, 1, 1);
  x_est_main.velocity = x_update.block(4, 0, 3, 1);
  x_est_main.acceleration = x_update.block(7, 0, 3, 1);
  x_est_main.bias = x_update.block(10, 0, num_of_gnss_channel, 1);
  //Cd = x_update(10);
  x_est_target.position = x_update.block(single_dimension, 0, 3, 1);
  x_est_target.clock = x_update.block(single_dimension + 3, 0, 1, 1);
  x_est_target.velocity = x_update.block(single_dimension + 4, 0, 3, 1);
  x_est_target.acceleration = x_update.block(single_dimension + 7, 0, 3, 1);
  x_est_target.bias = x_update.bottomRows(num_of_gnss_channel);

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(state_dimension, state_dimension);
  tmp = (I - K*H);
  M = tmp*M*tmp.transpose() + K*R*K.transpose(); // 負になるのが存在しているが，これは何？
  // M = (I - K * H) * M; // Mの観測更新がうまくいっていない．落ちすぎている．
  // TODO: Nについて，収束している？ものは0に落とす．
  for (int i = 0; i < num_of_gnss_channel; ++i)
  {
    int N_offset = num_of_single_status + i;
    if (M(N_offset, N_offset) < 0) M(N_offset, N_offset) = pow(sigma_N_ini, 2.0); // sqrt(pow(M(N_offset, N_offset), 2.0));
    N_offset = single_dimension + num_of_single_status + i;
    if (M(N_offset, N_offset) < 0) M(N_offset, N_offset) = pow(sigma_N_ini, 2.0); // sqrt(pow(M(N_offset, N_offset), 2.0));
  }
  // 数値誤差で発散してしまったやつを処理しないと．
  int n_main = main_info_.now_observed_gnss_sat_id.size();
  int offset = num_of_single_status + n_main;
  int size = num_of_gnss_channel - n_main;
  M.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);
  int n_target = target_info_.now_observed_gnss_sat_id.size();
  offset = single_dimension + num_of_single_status + n_target;
  size = num_of_gnss_channel - n_target;
  M.block(offset, offset, size, size) = Eigen::MatrixXd::Zero(size, size);
  
  // others
  for (int i = 0; i < num_of_single_status; ++i)
  {
    int offset = i;
    if (M(offset, offset) < 1e-6) M(offset, offset) = sqrt(pow(M(offset, offset), 2.0));
    offset = single_dimension + i;
    if (M(offset, offset) < 1e-6) M(offset, offset) = sqrt(pow(M(offset, offset), 2.0));
  }

  if (!std::isfinite(x_est_main.position(0)))
  {
    std::cout << "inf or nan" << std::endl;
    abort();
  }

  // clear
  for(int i=0; i<2;++i)
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
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(state_dimension, observation_dimension);
  K.topLeftCorner(num_of_single_status + observe_gnss_num_m, observe_gnss_num_m) = K_calc.topLeftCorner(num_of_single_status + observe_gnss_num_m, observe_gnss_num_m);
  K.block(single_dimension, 0, num_of_single_status + observe_gnss_num_t, observe_gnss_num_m) = K_calc.block(num_of_single_status + observe_gnss_num_m, 0, num_of_single_status + observe_gnss_num_t, observe_gnss_num_m);
  K.block(0, num_of_gnss_channel, num_of_single_status + observe_gnss_num_m, observe_gnss_num_t) = K_calc.block(0, observe_gnss_num_m, num_of_single_status + observe_gnss_num_m, observe_gnss_num_t);
  K.block(single_dimension, num_of_gnss_channel, num_of_single_status + observe_gnss_num_t, observe_gnss_num_t) = K_calc.block(num_of_single_status + observe_gnss_num_m, observe_gnss_num_m, num_of_single_status + observe_gnss_num_t, observe_gnss_num_t);
  K.block(0, 2*num_of_gnss_channel, num_of_single_status + observe_gnss_num_m, observe_gnss_num_c) = K_calc.block(0, observe_gnss_num_m+observe_gnss_num_t, num_of_single_status + observe_gnss_num_m, observe_gnss_num_c);
  K.block(single_dimension, 2*num_of_gnss_channel, num_of_single_status + observe_gnss_num_t, observe_gnss_num_c) = K_calc.block(num_of_single_status + observe_gnss_num_m, observe_gnss_num_m+observe_gnss_num_t, num_of_single_status + observe_gnss_num_t, observe_gnss_num_c);
  return K;
};

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
  RemoveRows(S, 2*num_of_gnss_channel + observe_gnss_c, 3*num_of_gnss_channel - 1);
  RemoveColumns(S, 2*num_of_gnss_channel + observe_gnss_c, 3*num_of_gnss_channel - 1);
  RemoveRows(S, num_of_gnss_channel + observe_gnss_t, 2*num_of_gnss_channel - 1);
  RemoveColumns(S, num_of_gnss_channel + observe_gnss_t, 2*num_of_gnss_channel - 1);
  RemoveRows(S, observe_gnss_m, num_of_gnss_channel - 1);
  RemoveColumns(S, observe_gnss_m, num_of_gnss_channel - 1);
};

void PBD_dgps::ResizeMHt(Eigen::MatrixXd& MHt, const int observe_gnss_m, const int observe_gnss_t, const int observe_gnss_c)
{
  // ここも後ろから．
  RemoveRows(MHt, single_dimension + num_of_single_status + observe_gnss_t, state_dimension - 1);
  RemoveRows(MHt, num_of_single_status + observe_gnss_m, single_dimension - 1);
  RemoveColumns(MHt, 2*num_of_gnss_channel + observe_gnss_c, 3*num_of_gnss_channel - 1);
  RemoveColumns(MHt, num_of_gnss_channel + observe_gnss_t, 2*num_of_gnss_channel - 1);
  RemoveColumns(MHt, observe_gnss_m, num_of_gnss_channel - 1);
};

void PBD_dgps::UpdateObservationsGRAPHIC(const int sat_id, EstimatedVariables& x_est, const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
{
  // ここもLEO satが把握している誤差ありの情報．
  // std::find index of the observing gnss satellite
  const GnssObserveInfo& observe_info_ = gnss_observations_.at(sat_id).info_;
  const int index = conv_index_from_gnss_sat_id(observe_info_.now_observed_gnss_sat_id, gnss_sat_id);
  const GnssObservedValues& observed_val_ = gnss_observations_.at(sat_id).observed_values_;
  auto gnss_position = observed_val_.gnss_satellites_position.at(index);
  double gnss_clock = observed_val_.gnss_clock.at(index);
  // とりあえずL1を使う．
  double pseudo_range = observed_val_.L1_pseudo_range.at(index);
  const auto& carrier_phase = observed_val_.L1_carrier_phase.at(index);
  double carrier_phase_range = (carrier_phase.first + carrier_phase.second) * L1_lambda;

  // ここら辺はGnssObserveModel内に格納する．
  double geometric_range = CalculateGeometricRange(x_est.position, gnss_position);
  gnss_observed_models_.at(sat_id).geometric_range.push_back(geometric_range);
  double pseudo_range_model = CalculatePseudoRange(x_est, gnss_position, gnss_clock);
  gnss_observed_models_.at(sat_id).pseudo_range_model.push_back(pseudo_range_model);
  double carrier_phase_range_model = CalculateCarrierPhase(x_est, gnss_position, gnss_clock, x_est.bias(index), L1_lambda);
  gnss_observed_models_.at(sat_id).carrier_phase_range_model.push_back(carrier_phase_range_model);

  // GRAPHIC
  const int row_offset = sat_id*num_of_gnss_channel + index;
  const int col_offset = sat_id*single_dimension + num_of_single_status + index;
  z(row_offset) = (pseudo_range + carrier_phase_range) / 2;
  h_x(row_offset) = (pseudo_range_model + carrier_phase_range_model) / 2;
  for (int j = 0; j < 3; ++j) {
    const int col_index = sat_id * single_dimension + j;
    // position
    H(row_offset, col_index) = (x_est.position(j) - gnss_position[j]) / geometric_range;
  }
  // clock
  H(row_offset, sat_id * single_dimension + 3) = 1.0;
  // Ambiguity
  H(row_offset, col_offset) = 0.5; // GRAPHIC N
  Rv(row_offset) = pow((pseudo_sigma + carrier_sigma) / 2.0, 2.0);
  // Rv(row_offset) = pow(pseudo_sigma, 2.0); // 少し大きめ
};

void PBD_dgps::UpdateObservationsSDCP(const int gnss_sat_id, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
{
  // std::find for comon index
  const int common_index = conv_index_from_gnss_sat_id(common_observed_gnss_sat_id, gnss_sat_id);
  const int row_offset = 2 * num_of_gnss_channel + common_index;

  const int main_index = conv_index_from_gnss_sat_id(gnss_observations_.at(0).info_.now_observed_gnss_sat_id, gnss_sat_id);
  // const int main_col_offset = 2 * num_of_gnss_channel + common_index;
  const int target_index = conv_index_from_gnss_sat_id(gnss_observations_.at(1).info_.now_observed_gnss_sat_id, gnss_sat_id);
  const int col_offset_main = num_of_single_status + main_index;
  const int col_offset_target = single_dimension + num_of_single_status + target_index;

  // とりあえずL1を使う．
  // main
  const auto& carrier_phase_main = gnss_observations_.at(0).observed_values_.L1_carrier_phase.at(main_index);
  double carrier_phase_range_main = (carrier_phase_main.first + carrier_phase_main.second) * L1_lambda;
  // target
  const auto& carrier_phase_target = gnss_observations_.at(1).observed_values_.L1_carrier_phase.at(target_index);
  double carrier_phase_range_target = (carrier_phase_target.first + carrier_phase_target.second) * L1_lambda;

  auto gnss_position = gnss_observations_.at(0).observed_values_.gnss_satellites_position.at(main_index);

  // SDCP
  z(row_offset) = carrier_phase_range_target - carrier_phase_range_main;
  h_x(row_offset) = gnss_observed_models_.at(1).carrier_phase_range_model.at(target_index) - gnss_observed_models_.at(0).carrier_phase_range_model.at(main_index);
  for (int j = 0; j < 3; ++j) {
    // position
    H(row_offset, j) = -(x_est_main.position(j) - gnss_position[j]) / gnss_observed_models_.at(0).geometric_range.at(main_index);
    H(row_offset, single_dimension + j) = (x_est_target.position(j) - gnss_position[j]) / gnss_observed_models_.at(1).geometric_range.at(target_index);
  }
  // clock
  H(row_offset, 3) = -1.0;
  H(row_offset, single_dimension + 3) = 1.0;
  // Ambiguity
  H(row_offset, col_offset_main) = -1.0; // SDCP N
  H(row_offset, col_offset_target) = 1.0; // SDCP N
  Rv(row_offset) = pow(sqrt(2.0) * carrier_sigma, 2.0);
};

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
void PBD_dgps::UpdateTrueBias(std::vector<std::vector<double>> bias, const int gnss_sat_id, const double lambda)
{
  const int index_main = conv_index_from_gnss_sat_id(gnss_observations_.at(0).info_.now_observed_gnss_sat_id, gnss_sat_id);
  true_bias_main(index_main) = bias[0].at(gnss_sat_id) * lambda;
  const int index_target = conv_index_from_gnss_sat_id(gnss_observations_.at(1).info_.now_observed_gnss_sat_id, gnss_sat_id);
  true_bias_target(index_target) = bias[1].at(gnss_sat_id) * lambda;
  /*
  if (now_main_observing_ch.count(gnss_sat_id))
  {
    int ch = now_main_observing_ch[gnss_sat_id];
    true_bias_main(ch) = bias[0].at(gnss_sat_id)*lambda;
  }
  if (now_common_observing_ch.count(gnss_sat_id))
  {
    int ch = now_common_observing_ch[gnss_sat_id];
    true_bias_target(ch) = bias[1].at(gnss_sat_id)*lambda;
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

  const Eigen::VectorXd pre_estimated_bias = x_est.bias; // 参照渡し？
  Eigen::MatrixXd pre_M = M;
  // reset 
  M.block(sat_id*single_dimension + num_of_single_status, sat_id*single_dimension + num_of_single_status, num_of_gnss_channel, num_of_gnss_channel) = Eigen::MatrixXd::Zero(num_of_gnss_channel, num_of_gnss_channel);

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
      x_est.bias(now_index) = N_dist(mt); // こいつもGauss Makov過程にした方がいいらしい?
      int offset = sat_id * single_dimension + num_of_single_status + now_index;
      // 行と列に関してもリセット
      M.block(0, offset, state_dimension, 1) = Eigen::MatrixXd::Zero(state_dimension, 1);
      M.block(offset, 0, 1, state_dimension) = Eigen::MatrixXd::Zero(1, state_dimension);
      M(offset, offset) = pow(sigma_N_ini, 2.0);
      ++now_index;
    }
    // 引き継ぐ
    else if (observe_info_.pre_observed_status.at(i) == true && observe_info_.now_observed_status.at(i) == true)
    {
      x_est.bias(now_index) = pre_estimated_bias(pre_index);
      int offset_base = sat_id * single_dimension + num_of_single_status;
      M.block(0, offset_base + now_index, state_dimension, 1) = pre_M.block(0, offset_base + pre_index, state_dimension, 1);
      M.block(offset_base + now_index, 0, 1, state_dimension) = pre_M.block(offset_base + pre_index, 0, 1, state_dimension);
      ++pre_index; ++now_index;
    }
    if (now_index >= num_of_gnss_channel || pre_index >= num_of_gnss_channel) break; // ch以上の受信は出来ない
  }
  // now_index以上の部分を0に落とすということをやる．
  x_est.bias.block(now_index, 0, num_of_gnss_channel - now_index, 1) = Eigen::VectorXd::Zero(num_of_gnss_channel - now_index);


  /*
  if (abs(M.determinant()) < 10e-13)
  {
    cout << "M matirx is singular" << std::endl;
    abort();
  }
  */
  //cout << M << std::endl;
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


double PBD_dgps::CalculatePseudoRange(const EstimatedVariables& x_est, libra::Vector<3> gnss_position, double gnss_clock) // 本来はここのgnss_positionをIGSとかとIGUとかで分けて扱う必要がある
{
    double res = 0.0;
    res = CalculateGeometricRange(x_est.position, gnss_position);
    // clock offsetの分を追加
    res += x_est.clock(0) - gnss_clock; // 電離層はフリーにしている．
    // 観測ノイズ
    std::normal_distribution<> pseudo_range_noise(0.0, pseudo_sigma);
    res += pseudo_range_noise(mt);
    return res;
}

// 電離圏とかはまた考える．<- ここで追加できないなら電離圏フリーの観測量じゃないとダメ．
double PBD_dgps::CalculateCarrierPhase(const EstimatedVariables& x_est, libra::Vector<3> gnss_position, double gnss_clock, double integer_bias, double lambda)
{
  double res = 0.0;
  res = CalculateGeometricRange(x_est.position, gnss_position);
  res += x_est.clock(0) - gnss_clock;
  // res -= lambda * integer_bias; // この形でできるようにするためにはbiasがちゃんと整数不定性の量になっている必要がある．
  res += integer_bias;
  // 観測ノイズ
  std::normal_distribution<> carrier_phase_noise(0.0, carrier_sigma);
  res += carrier_phase_noise(mt);

  return res;
}

double PBD_dgps::CalculateGeometricRange(const Eigen::Vector3d& sat_position, libra::Vector<3> gnss_position) const
{
    double res = 0.0;
    for(int i = 0;i < 3;++i){
        res += pow(sat_position(i) - gnss_position[i], 2.0);
    }
    res = sqrt(res);

    return res;
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