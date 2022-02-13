#include <iomanip>
#include <set>
#include "PBD_dgps.h"
#include "PBD_const.h"

static const int conv_index_from_gnss_sat_id(vector<int> observed_gnss_sat_id, const int gnss_sat_id);

// outputを変えるときは"result.csv"を変更する．せめてパスは変えたい．
PBD_dgps::PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Orbit& main_orbit, const Orbit& target_orbit):mt(42), sat_main_clock_true(0.0), sat_target_clock_true(0.0), step_time(sim_time_.GetStepSec()), ofs("result_new.csv"), num_of_gnss_satellites(gnss_satellites_.GetNumOfSatellites()), main_orbit_(main_orbit), target_orbit_(target_orbit)
{
  //初期化
  x_est_main.position = Eigen::VectorXd::Zero(3);
  libra::Vector<3> position_main = main_orbit_.GetSatPosition_i();
  for(int i = 0;i < 3;++i) x_est_main.position(i) = position_main[i];
  x_est_main.clock = Eigen::VectorXd::Zero(1);
  libra::Vector<3> velocity_main = main_orbit_.GetSatVelocity_i();
  for(int i = 0;i < 3;++i) x_est_main.velocity(i) = velocity_main[i];
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

  // x_est = { x_est_main, x_est_target };
  // ここはpush_backで問題なかったっポイ
  //x_est.reserve(2);
  x_est.push_back(x_est_main);
  x_est.push_back(x_est_target);
  //x_est.at(0) = x_est_main;
  //x_est.at(1) = x_est_target;
  // observe_info.reserve(2);
  // ここが参照になっていないのがつらい．
  observe_info.push_back(observe_info_main);
  observe_info.push_back(observe_info_target);
  // observe_info.at(0) = observe_info_main;
  // observe_info.at(1) = observe_info_target;

  true_bias_main = Eigen::VectorXd::Zero(num_of_gnss_channel);
  true_bias_target = Eigen::VectorXd::Zero(num_of_gnss_channel);

  // 初期分散
  std::normal_distribution<> position_dist(0.0,sigma_r_ini);
  std::normal_distribution<> clock_dist(0.0, sigma_cdt_ini);
  std::normal_distribution<> velocity_dist(0.0, sigma_v_ini);
  std::normal_distribution<> acc_dist(0.0, sigma_acc_r_ini);
  std::normal_distribution<> N_dist(0.0, sigma_N_ini);

  // Vって名前は変えた方がいいかも
  Eigen::VectorXd V = Eigen::VectorXd::Constant(state_dimension, 0);

  // A-priori 
  for(int i = 0; i < 3; ++i) V(i) = pow(sigma_r_ini, 2.0); // position
  V(3) = pow(sigma_cdt_ini, 2.0); // clock
  for(int i = 0; i < 3; ++i) V(4 + i) = pow(sigma_v_ini, 2.0); // velocity
  for(int i = 0; i < 3; ++i) V(7 + i) = pow(sigma_acc_r_ini, 2.0); // acceleration RTNになっていないので後ほど修正
  // 初期は全部に入れている．
  for (int i = 0; i < num_of_gnss_channel; ++i) V(num_of_single_status + i) = pow(sigma_N_ini, 2.0); // N
  // 以下がtarget
  for (int i = 0; i < 3; ++i) V(single_dimension + i) = pow(sigma_r_ini, 2.0);
  V(single_dimension + 3) = pow(sigma_cdt_ini, 2.0);
  for (int i = 0; i < 3; ++i) V(single_dimension + 4 + i) = pow(sigma_v_ini, 2.0);
  for (int i = 0; i < 3; ++i) V(single_dimension + 7 + i) = pow(sigma_acc_r_ini, 2.0); // ここもRTNになってないので後ほど修正
  for (int i = 0; i < num_of_gnss_channel; ++i) V(single_dimension + num_of_single_status + i) = pow(sigma_N_ini, 2.0); // N

  M = V.asDiagonal(); // 誤差共分散行列M

  // 初期位置はガウシアンからサンプリング．mtは乱数のシード
  for(int i = 0; i < 3; ++i) x_est_main.position(i) += position_dist(mt);
  x_est_main.clock(0) += clock_dist(mt);
  for(int i = 0; i < 3; ++i) x_est_main.velocity(i) += velocity_dist(mt);
  for(int i = 0; i < 3; ++i) x_est_main.acceleration(i) += acc_dist(mt);
  for(int i = 0; i < num_of_gnss_channel; ++i) x_est_main.bias(i) += N_dist(mt);
  for(int i = 0; i < 3; ++i) x_est_target.position(i) += position_dist(mt);
  x_est_target.clock(0) += clock_dist(mt);
  for (int i = 0; i < 3; ++i) x_est_target.velocity(i) += velocity_dist(mt);
  for(int i = 0; i < 3; ++i) x_est_target.acceleration(i) += acc_dist(mt);
  for(int i = 0; i < num_of_gnss_channel; ++i) x_est_target.bias(i) += N_dist(mt);

  for (int i = 0; i < 2; ++i)
  {
    observe_info.at(i).pre_observed_status.assign(num_of_gnss_satellites, false);
    observe_info.at(i).now_observed_status.assign(num_of_gnss_satellites, false);
    L1_bias.at(i).assign(num_of_gnss_satellites, 0.0);
    L2_bias.at(i).assign(num_of_gnss_satellites, 0.0);
  }
  common_observed_status.assign(num_of_gnss_satellites, false);
  for (int i = 0; i < num_of_gnss_channel; ++i) main_free_ch.push_back(i);
  for (int i = 0; i < num_of_gnss_channel; ++i) common_free_ch.push_back(i);

  ofstream ofs_ini_txt("readme_new.txt");
  ofs_ini_txt << "initial position dist: " << sigma_r_ini << endl;
  ofs_ini_txt << "initial velocity dist: " << sigma_v_ini << endl;
  ofs_ini_txt << "initial acceleration dist: " << sigma_acc_r_ini << endl;
  ofs_ini_txt << "initial clock dist: " << sigma_cdt_ini << endl;
  ofs_ini_txt << "initial ambiguity dist: " << sigma_N_ini << endl;
  ofs_ini_txt << "pseudo dist: " << pseudo_sigma << endl;
  ofs_ini_txt << "carrier dist: " << carrier_sigma << endl;
  ofs_ini_txt << "clock dist: " << clock_sigma << endl;
  ofs_ini_txt << "process noise of acceleration: " << sigma_acc_r_process << endl;
  ofs_ini_txt << "process noise of clock: " << sigma_cdt_process << endl;
  ofs_ini_txt << "process noise of ambiguity: " << sigma_N_process << endl;
  ofs_ini_txt << "mask angle: " << mask_angle << endl;
  ofs_ini_txt << "num of status: " << num_of_status << endl;
  ofs_ini_txt << "observe step time: " << observe_step_time << endl;
  ofs_ini_txt << "log step time: " << log_step_time << endl;
}

PBD_dgps::~PBD_dgps(){}

void PBD_dgps::Update(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_)// , const Orbit& main_orbit, const Orbit& target_orbit)
{
  //clock //仮
  std::normal_distribution<> clock_dist(0.0, clock_sigma);
  //  この扱いが意味わからん
  sat_main_clock_true = clock_dist(mt); // main
  sat_target_clock_true = clock_dist(mt); // target

  double elapsed_time = sim_time_.GetElapsedSec();
  double tmp = floor(elapsed_time/observe_step_time + 1e-4); //1e-4は数値誤差
  double tmp_log = floor(elapsed_time/log_step_time + 1e-4);

  //まず更新
  OrbitPropagation();

  if(abs(elapsed_time - tmp*observe_step_time) < step_time/2.0){
    //観測時間にピッタリ
    GnssObservedValues gnss_observed_values_main;
    GnssObservedValues gnss_true_main;
    GetGnssPositionObservation(gnss_satellites_, main_orbit_, 0, gnss_observed_values_main, gnss_true_main, sat_main_clock_true);
    // Ionfreeを生成している．
    ProcessGnssObservation(gnss_observed_values_main, gnss_true_main, 0);
    GnssObservedValues gnss_observed_values_target;
    GnssObservedValues gnss_true_target;
    GetGnssPositionObservation(gnss_satellites_, target_orbit_, 1, gnss_observed_values_target, gnss_true_target, sat_target_clock_true);
    ProcessGnssObservation(gnss_observed_values_target, gnss_true_target, 1);
    //std::pair<GnssObservedValues, GnssObservedValues> gnss_observed_pair = std::make_pair(gnss_observed_values_main, gnss_observed_values_target);
    //std::pair<GnssObservedValues, GnssObservedValues> gnss_true_pair = std::make_pair(gnss_true_main, gnss_true_target);

    // ここにFindCommonSatellite?
    SetBiasToObservation(0, x_est_main, gnss_observed_values_main, gnss_true_main);
    SetBiasToObservation(1, x_est_target, gnss_observed_values_target, gnss_true_target);
    KalmanFilter(gnss_observed_values_main, gnss_observed_values_target);
  }

  //log output
  if(abs(elapsed_time - tmp_log*log_step_time) < step_time/2.0){
    libra::Vector<3> sat_position = main_orbit_.GetSatPosition_i();
    libra::Vector<3> sat_velocity = main_orbit_.GetSatVelocity_i();
    for(int i = 0;i < 3;++i) ofs << fixed << setprecision(30) << sat_position[i] << ","; // r_m_true
    ofs << fixed << setprecision(30) << sat_main_clock_true << ","; // t_m_true
    for(int i = 0;i < 3;++i) ofs << fixed << setprecision(30) << sat_velocity[i] << ","; // v_m_ture
    libra::Vector<3> sat_position_target = target_orbit_.GetSatPosition_i();
    libra::Vector<3> sat_velocity_target = target_orbit_.GetSatVelocity_i();
    for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << sat_position_target[i] << ","; // r_t_ture
    ofs << fixed << setprecision(30) << sat_target_clock_true << ","; // t_t_true
    for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << sat_velocity_target[i] << ","; // v_t_ture
    for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << x_est_main.position(i) << ","; // r_m_est
    ofs << fixed << setprecision(30) << x_est_main.clock(0) << ","; // t_m_est
    for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << x_est_main.velocity(i) << ","; // v_m_est
    for(int i = 0;i < 3;++i) ofs << fixed << setprecision(30) << x_est_main.acceleration(i) << ","; // a_m_est
    for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << x_est_target.position(i) << ","; // r_t_est
    ofs << fixed << setprecision(30) << x_est_target.clock(0) << ","; //t_t_est
    for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << x_est_target.velocity(i) << ","; // v_t_est
    for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << x_est_target.acceleration(i) << ","; // a_t_est
    for (int i = 0; i < num_of_gnss_channel; ++i) ofs << fixed << setprecision(30) << true_bias_main(i) << ","; // N_true
    for (int i = 0; i < num_of_gnss_channel; ++i) ofs << fixed << setprecision(30) << x_est_main.bias(i) << ","; // N_est
    for (int i = 0; i < num_of_gnss_channel; ++i) ofs << fixed << setprecision(30) << true_bias_target(i) << ","; // dN_true
    for (int i = 0; i < num_of_gnss_channel; ++i) ofs << fixed << setprecision(30) << x_est_target.bias(i) << ","; // dN_est
    for (int i = 0; i < state_dimension; ++i) ofs << fixed << setprecision(30) << M(i, i) << ",";
    // ここは何を残しているのか確認．
    int ans = 0;
    for(int i = 0;i < num_of_gnss_satellites;++i){
        auto gnss_position = gnss_satellites_.GetSatellitePositionEci(i);
        if(CheckCanSeeSatellite(sat_position, gnss_position)) ++ans;
    }
    ofs << ans << endl;
  }

  return;
}

// aがおかしい気がする．単位の間違いか？
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
  vector<Eigen::Vector3d> pos_vel_main = RK4(position_main, velocity_main, acceleration_main);
  position_main = pos_vel_main[0];
  velocity_main = pos_vel_main[1];

  vector<Eigen::Vector3d> pos_vel_target = RK4(position_target, velocity_target, acceleration_target);
  position_target = pos_vel_target[0];
  velocity_target = pos_vel_target[1];

  x_est_main.position = position_main;
  x_est_main.velocity = velocity_main;
  x_est_target.position = position_target;
  x_est_target.velocity = velocity_target;

  // acceleration ここでは q = sigma_acc_process
  Eigen::MatrixXd Phi_a = CalculatePhi_a(step_time);
  double phi = Phi_a(0, 0);
  x_est_main.acceleration = Phi_a * x_est_main.acceleration;
  x_est_target.acceleration = Phi_a * x_est_target.acceleration;
  std::normal_distribution<> acc_r_process_noise(0.0, sigma_acc_r_process*sqrt(1 - phi*phi));
  for (int i = 0; i < 3; ++i) x_est_main.acceleration(i) += acc_r_process_noise(mt);
  for (int i = 0; i < 3; ++i) x_est_target.acceleration(i) += acc_r_process_noise(mt);
  // cdt
  std::normal_distribution<> cdt_process_noise(0.0, sigma_cdt_process*sqrt(step_time/tau_cdt));
  x_est_main.clock(0) += cdt_process_noise(mt);
  x_est_target.clock(0) += cdt_process_noise(mt);

  M = UpdateM();
}

// Orbitの更新を使えるように修正．位置と速度はRK4で伝搬しているのが，共分散はオイラー法になっている．EKFなので仕方がない．
vector<Eigen::Vector3d> PBD_dgps::RK4(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration)
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

  vector<Eigen::Vector3d> position_velocity = 
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

  double ac_norm = -mu_const/position.squaredNorm(); //2体の重力項
  double tmp_J2_coefficient = 3.0/2.0*mu_const*J2_const*pow(Earth_Radius, 2.0)/pow(r, 4.0); //J2項の係数

  Eigen::Vector3d all_acceleration = position/r;

  all_acceleration(0) *= ac_norm - tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0));
  all_acceleration(1) *= ac_norm - tmp_J2_coefficient*(1.0 - 5.0*pow(z/r, 2.0));
  all_acceleration(2) *= ac_norm - tmp_J2_coefficient*(3.0 - 5.0*pow(z/r, 2.0));

  all_acceleration -= Cd*v*velocity; //-Cd*V^2*(Vi/V) 大気抵抗

  all_acceleration += acceleration*conv_nm2m; //残りの摂動要素

  return all_acceleration; // m/s2
}

// Mの各要素の単位を改善した方がいい．急激に小さい量だけに減るのでMの正則性が失われている．
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
  // clock
  Phi(3, 3) = 1.0; // random walkなのでΦは1
  Phi(single_dimension + 3, single_dimension + 3) = 1.0;
  Phi += step_time * A;
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
  Eigen::MatrixXd Gamma = pow(step_time, 2.0) * B * Q * B.transpose(); //pow(step_time, 2.0)?

  // 他のprocess noiseも入れて感度を見てみる．
  Eigen::MatrixXd res = Phi * M * Phi.transpose() + Gamma;
  // Nのない部分を0に落とす？
  int n_main = observe_info.at(0).now_observed_gnss_sat_id.size();
  res.block(num_of_single_status + n_main, num_of_single_status + n_main, num_of_gnss_channel - n_main, num_of_gnss_channel - n_main) = Eigen::MatrixXd::Zero(num_of_gnss_channel - n_main, num_of_gnss_channel - n_main);
  int n_target = observe_info.at(1).now_observed_gnss_sat_id.size();
  res.block(single_dimension + num_of_single_status + n_target, single_dimension + num_of_single_status + n_target, num_of_gnss_channel - n_target, num_of_gnss_channel - n_target) = Eigen::MatrixXd::Zero(num_of_gnss_channel - n_target, num_of_gnss_channel - n_target);
  
  return res;
}

/*
// Jacobi or STM for pod
Eigen::MatrixXd PBD_dgps::CalculateA(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration) const
{
    double r = position.norm();
    double v = velocity.norm();

    double x = position(0); double y = position(1); double z = position(2);
    double vx = velocity(0); double vy = velocity(1); double vz = velocity(2);

    double J2_coefficient = 3.0/2.0*mu_const*J2_const*Earth_Radius*Earth_Radius; 

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_of_status, num_of_status);
    // (r, v)
    A(0,4) = 1.0; A(1,5) = 1.0; A(2,6) = 1.0;

    // (v, r)
    A(4,0) = 3.0*mu_const*x*x/pow(r, 5.0) - mu_const/pow(r, 3.0) - J2_coefficient*(1.0/pow(r, 5.0) - 5.0*(x*x + z*z)/pow(r, 7.0) + 35.0*x*x*z*z/pow(r, 9.0));
    A(4,1) = 3.0*mu_const*x*y/pow(r, 5.0) - J2_coefficient*( -5.0*x*y/pow(r, 7.0) + 35.0*x*y*z*z/pow(r, 9.0));
    A(4,2) = 3.0*mu_const*x*z/pow(r, 5.0) - J2_coefficient*( -15.0*x*z/pow(r, 7.0) + 35.0*x*z*z*z/pow(r, 9.0));
    A(5,0) = 3.0*mu_const*x*y/pow(r, 5.0) - J2_coefficient*( -5.0*x*y/pow(r, 7.0) + 35.0*x*y*z*z/pow(r, 9.0));
    A(5,1) = 3.0*mu_const*y*y/pow(r, 5.0) - mu_const/pow(r, 3.0) - J2_coefficient*(1.0/pow(r, 5.0) - 5.0*(y*y + z*z)/pow(r, 7.0) + 35.0*y*y*z*z/pow(r, 9.0));
    A(5,2) = 3.0*mu_const*y*z/pow(r, 5.0) - J2_coefficient*( -15.0*y*z/pow(r, 7.0) + 35.0*y*z*z*z/pow(r, 9.0));
    A(6,0) = 3.0*mu_const*x*z/pow(r, 5.0) - J2_coefficient*( -15.0*x*z/pow(r, 7.0) + 35.0*x*z*z*z/pow(r, 9.0));
    A(6,1) = 3.0*mu_const*y*z/pow(r, 5.0) - J2_coefficient*( -15.0*y*z/pow(r, 7.0) + 35.0*y*z*z*z/pow(r, 9.0));
    A(6,2) = 3.0*mu_const*z*z/pow(r, 5.0) - mu_const/pow(r, 3.0) - J2_coefficient*(3.0/pow(r, 5.0) - 30.0*z*z/pow(r, 7.0) + 35.0*pow(z, 4.0)/pow(r, 9.0));

    // これ速度に入る？普通は加速度として入ってくるはず．
    // A(4,4) = -Cd*(vx*vx/v + v);    A(4,5) = -Cd*vx*vy/v;    A(4,6) = -Cd*vx*vz/v;
    // A(5,4) = -Cd*vx*vy/v;    A(5,5) = -Cd*(vy*vy/v + v);    A(5,6) = -Cd*vy*vz/v;
    // A(6,4) = -Cd*vx*vz/v;    A(6,5) = -Cd*vy*vz/v;    A(6,6) = -Cd*(vz*vz/v + v);
    
    // これを入れるなら1/2t^2も入れたい<-非線形こうなので無理？
    // A(4,7) = 1.0;	A(5,8) = 1.0;	A(6,9) = 1.0;

    //A(4,10) = -v*vx;    A(5,10) = -v*vy;    A(6,10) = -v*vz;

    return A;
}
*/

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
  A(4, 7) = 1.0*1e-6;	A(5, 8) = 1.0 * 1e-6;	A(6, 9) = 1.0 * 1e-6;

  // これはもしCdも推定しているならいる．
  //A(4,10) = -v*vx;    A(5,10) = -v*vy;    A(6,10) = -v*vz;
  return A;
};

// Process noiseのvarianceを計算．Bを使う形に修正．
Eigen::MatrixXd PBD_dgps::CalculateQ()
{
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(8, 8);
  double phi = exp(-step_time / tau_a);
  double q_acc = pow(sigma_acc_r_process, 2.0) * (1 - pow(phi, 2.0));
  double q_cdt = sigma_cdt_process * sqrt(step_time / tau_cdt);
  // for(int i = 0;i < 3;++i) Q(i, i) = pow(sigma_r_process, 2.0);
  Q(0, 0) = q_cdt;
  // for(int i = 4;i < 7;++i) Q(i, i) = pow(sigma_v_process, 2.0);
  for(int i = 0;i < 3;++i) Q(1 + i, 1 + i) = q_acc;
  // for (int i = num_of_single_status; i < n_main + num_of_single_status; ++i) Q(i, i) = pow(sigma_N_process, 2.0);
  // for (int i = num_of_single_status; i < num_of_single_status + 3; ++i) Q(i, i) = pow(sigma_r_process, 2.0);
  Q(4, 4) = q_cdt;
  // for (int i = num_of_single_status + 4; i < num_of_single_status + 7; ++i) Q(i, i) = pow(sigma_v_process, 2.0);
  for(int i = 0; i < 3; ++i) Q(5 + i, 5 + i) = q_acc;
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

// 観測量と状態量をchで管理するのか上から埋めていくのかはよく考える．数学的に問題がないのかを
// 行列サイズ一定方針は無理なんかも．
void PBD_dgps::KalmanFilter(const GnssObservedValues& gnss_observed_main, const GnssObservedValues& gnss_observed_target)
{
  // gnss_observed_hoge の観測される時間はどうなっているのか？
  // int n_common = common_observed_gnss_sat_id.size();
  
  // GRAPHIC*2 + SDCPにする．
  Eigen::VectorXd z = Eigen::VectorXd::Zero(observation_dimension); //観測ベクトル
  Eigen::VectorXd h_x = Eigen::VectorXd::Zero(observation_dimension); // 観測モデル行列
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(observation_dimension, state_dimension); //観測行列（dhx/dx）
  // Oneにしている部分がどう効いてくるのかの確認は必要．
  Eigen::VectorXd R_V = Eigen::MatrixXd::Ones(observation_dimension, 1); // 観測誤差共分散．
  // 正則になるように加工
  // for (int i = 0; i < num_of_gnss_channel; ++i) H.block(i, 0, 1, 3) = Eigen::MatrixXd::Ones(1, 3);
  // for (int i = 0; i < num_of_gnss_channel; ++i) H.block(num_of_gnss_channel + i, 0, 1, 3) = Eigen::MatrixXd::Ones(1, 3);
  // ここの入れ方も汎用的になるように関数かしてしまうのがいいかも．ただ，main, target, commonでそれぞれ回すと計算の無駄が多い．main targetの観測している衛星全部をまとめたものがあればいいんか．
  vector<int> all_observed_gnss_ids = observe_info.at(0).now_observed_gnss_sat_id;
  all_observed_gnss_ids.insert(all_observed_gnss_ids.end(), observe_info.at(1).now_observed_gnss_sat_id.begin(), observe_info.at(1).now_observed_gnss_sat_id.end()); // concate
  sort(all_observed_gnss_ids.begin(), all_observed_gnss_ids.end());
  all_observed_gnss_ids.erase(unique(all_observed_gnss_ids.begin(), all_observed_gnss_ids.end()), all_observed_gnss_ids.end()); // unique
  // この時のcommonがあるものから並んでいるとうれしい．<- 順番に関してはいったんかんがえんくていいか．
  for (const int& id :all_observed_gnss_ids)
  {
    // if main
    // int id = all_observed_gnss_ids.at(i);
    if (*find(observe_info.at(0).now_observed_gnss_sat_id.begin(), observe_info.at(0).now_observed_gnss_sat_id.end(), id) == id)
    {
      UpdateObservationsGRAPHIC(0, x_est_main, id, gnss_observed_main, z, h_x, H, R_V);
    }
    // if target
    if (*find(observe_info.at(1).now_observed_gnss_sat_id.begin(), observe_info.at(1).now_observed_gnss_sat_id.end(), id) == id)
    {
      UpdateObservationsGRAPHIC(1, x_est_target, id, gnss_observed_target, z, h_x, H, R_V);
    }
    // if common
    if (*find(common_observed_gnss_sat_id.begin(), common_observed_gnss_sat_id.end(), id) == id)
    {
      UpdateObservationsSDCP(id, gnss_observed_main, gnss_observed_target, z, h_x, H, R_V);
    }
  }
  // R_V.block(0, 0, 2*num_of_gnss_channel, 1) = Eigen::VectorXd::Constant(2*num_of_gnss_channel, pseudo_sigma / 2.0);
  // R_V.block(2*num_of_gnss_channel, 0, num_of_gnss_channel, 1) = Eigen::VectorXd::Constant(num_of_gnss_channel, sqrt(2.0) * carrier_sigma);
  // clock
  // H.block(0, 3, num_of_gnss_channel, 1) = Eigen::MatrixXd::Ones(num_of_gnss_channel, 1);
    
  Eigen::MatrixXd R = R_V.asDiagonal();

  Eigen::MatrixXd hmh = H * M * H.transpose();
  /*
  if (abs(hmh.determinant()) < 10e-10)
  {
    cout << "HMHt matirx is singular" << endl;
    abort();
  }
  */
  Eigen::MatrixXd tmp = R + hmh; // (observation_num, observation_num)
  // 観測量がないところは0にし直すみたいな加工が必要かもしれない．数字の桁数で打ち切りみたいな形にできればいい．
  // 正則チェック
  if (abs(tmp.determinant()) < 1e-15)
  {
    cout << "tmp matirx is singular" << endl; //めちゃくちゃ引っかかっている．
    // abort();
  }
  //カルマンゲイン
  Eigen::MatrixXd K = M * H.transpose() * tmp.inverse();
    
  // 観測量の無い項の影響をなくすために
  /*
  int n_main = observe_info_main.now_observed_gnss_sat_id.size();
  for (int i=0; i<num_of_gnss_channel-n_main; ++i)
  {
    K.block(0, n_main + i, state_dimension, 1) = Eigen::MatrixXd::Zero(state_dimension, 1);
  }
  int n_target = observe_info_target.now_observed_gnss_sat_id.size();
  for (int i = 0; i < num_of_gnss_channel - n_target; ++i)
  {
    K.block(0, num_of_gnss_channel + n_target + i, state_dimension, 1) = Eigen::MatrixXd::Zero(state_dimension, 1);
  }
  int n_common = common_observed_gnss_sat_id.size();
  for (int i = 0; i < num_of_gnss_channel - n_common; ++i)
  {
    K.block(0, 2*num_of_gnss_channel + n_common + i, state_dimension, 1) = Eigen::MatrixXd::Zero(state_dimension, 1);
  }
  */

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

  M = M - K * H * M; // rに対するMの観測更新がうまくいっていない．
  //Eigen::MatrixXd H_trans = H.transpose();
  //M = M - H.inverse()*R*H_trans.inverse();

  if (!std::isfinite(x_est_main.position(0)))
  {
    cout << "inf or nan" << endl;
    abort();
  }
  // clear
  for(int i=0; i<2;++i)
  {
    observe_info.at(i).geometric_range.clear();
    observe_info.at(i).pseudo_range_model.clear();
    observe_info.at(i).carrier_phase_range_model.clear();
  }
  return;
}

void PBD_dgps::UpdateObservationsGRAPHIC(const int sat_id, EstimatedVariables& x_est, const int gnss_sat_id, const GnssObservedValues& gnss_observed_val, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
{
  // ここもLEO satが把握している誤差ありの情報．
  // find index of the observing gnss satellite
  const int index = conv_index_from_gnss_sat_id(observe_info.at(sat_id).now_observed_gnss_sat_id, gnss_sat_id);
  auto gnss_position = gnss_observed_val.gnss_satellites_position.at(index);
  double gnss_clock = gnss_observed_val.gnss_clock.at(index);
  // とりあえずL1を使う．
  double pseudo_range = gnss_observed_val.L1_pseudo_range.at(index);
  const auto& carrier_phase = gnss_observed_val.L1_carrier_phase.at(index);
  double carrier_phase_range = (carrier_phase.first + carrier_phase.second) * L1_lambda;


  double geometric_range = CalculateGeometricRange(x_est.position, gnss_position);
  observe_info.at(sat_id).geometric_range.push_back(geometric_range);
  double pseudo_range_model = CalculatePseudoRange(x_est, gnss_position, gnss_clock);
  observe_info.at(sat_id).pseudo_range_model.push_back(pseudo_range_model);
  double carrier_phase_range_model = CalculateCarrierPhase(x_est, gnss_position, gnss_clock, x_est.bias(index), L1_lambda);
  observe_info.at(sat_id).carrier_phase_range_model.push_back(carrier_phase_range_model);

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
};

void PBD_dgps::UpdateObservationsSDCP(const int gnss_sat_id, const GnssObservedValues& gnss_observed_main, const GnssObservedValues& gnss_observed_target, Eigen::VectorXd& z, Eigen::VectorXd& h_x, Eigen::MatrixXd& H, Eigen::VectorXd& Rv)
{
  // find for comon index
  const int common_index = conv_index_from_gnss_sat_id(common_observed_gnss_sat_id, gnss_sat_id);
  const int row_offset = 2 * num_of_gnss_channel + common_index;

  const int main_index = conv_index_from_gnss_sat_id(observe_info.at(0).now_observed_gnss_sat_id, gnss_sat_id);
  // const int main_col_offset = 2 * num_of_gnss_channel + common_index;
  const int target_index = conv_index_from_gnss_sat_id(observe_info.at(1).now_observed_gnss_sat_id, gnss_sat_id);
  const int col_offset_main = num_of_single_status + main_index;
  const int col_offset_target = single_dimension + num_of_single_status + target_index;

  // とりあえずL1を使う．
  // main
  const auto& carrier_phase_main = gnss_observed_main.L1_carrier_phase.at(main_index);
  double carrier_phase_range_main = (carrier_phase_main.first + carrier_phase_main.second) * L1_lambda;
  // target
  const auto& carrier_phase_target = gnss_observed_target.L1_carrier_phase.at(target_index);
  double carrier_phase_range_target = (carrier_phase_target.first + carrier_phase_target.second) * L1_lambda;
  
  auto gnss_position = gnss_observed_main.gnss_satellites_position.at(main_index);

  // SDCP
  z(row_offset) = carrier_phase_range_target - carrier_phase_range_main;
  h_x(row_offset) = observe_info.at(1).carrier_phase_range_model.at(target_index) - observe_info.at(0).carrier_phase_range_model.at(main_index);
  for (int j = 0; j < 3; ++j) {
    // position
    H(row_offset, j) = -(x_est_main.position(j) - gnss_position[j]) / observe_info.at(0).geometric_range.at(main_index);
    H(row_offset, single_dimension + j) = (x_est_main.position(j) - gnss_position[j]) / observe_info.at(1).geometric_range.at(target_index);
  }
  // clock
  H(row_offset, 3) = -1.0;
  H(row_offset, single_dimension + 3) = 1.0;
  // Ambiguity
  H(row_offset, col_offset_main) = -1.0; // SDCP N
  H(row_offset, col_offset_target) = 1.0; // SDCP N
  Rv(row_offset) = pow(sqrt(2.0)*carrier_sigma, 2.0);
};

static const int conv_index_from_gnss_sat_id(vector<int> observed_gnss_sat_id, const int gnss_sat_id)
{
  vector<int>::iterator itr = find(observed_gnss_sat_id.begin(), observed_gnss_sat_id.end(), gnss_sat_id);
  if (itr == observed_gnss_sat_id.end())
  {
    cout << "not found" << gnss_sat_id << endl;
    abort();
  }
  const int index = distance(observed_gnss_sat_id.begin(), itr);
  return index;
};

void PBD_dgps::GetGnssPositionObservation(const GnssSatellites& gnss_satellites_, const Orbit& orbit_, const int sat_id, GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, double sat_clock_true)
{
    //推定値の計算
    int num_of_gnss_satellites = gnss_satellites_.GetNumOfSatellites();
    libra::Vector<3> sat_position_i = orbit_.GetSatPosition_i(); // ECI?

    observe_info.at(sat_id).now_observed_gnss_sat_id.clear(); //クリア

    for(int i = 0; i < num_of_gnss_satellites; ++i){
        //if(i == 7 || i == 23 || i == 31) continue;
        if(!gnss_satellites_.GetWhetherValid(i)) continue;
        libra::Vector<3> gnss_position = gnss_satellites_.Get_true_info().GetSatellitePositionEci(i);
        bool see_flag = CheckCanSeeSatellite(sat_position_i, gnss_position);
        // init
        main_index_dict.insert(std::make_pair(i, -1));
        common_index_dict.insert(std::make_pair(i, -1));

        int gnss_sat_id = i;
        if (!see_flag)
        {
          pre_main_observing_ch = now_main_observing_ch;
          RemoveFromCh(gnss_sat_id, now_main_observing_ch, main_free_ch);
          continue;
        }
        observe_info.at(sat_id).now_observed_status.at(gnss_sat_id) = true;
        observe_info.at(sat_id).now_observed_gnss_sat_id.push_back(gnss_sat_id);
        int observed_gnss_index = observe_info.at(sat_id).now_observed_gnss_sat_id.size() - 1;
        if (sat_id == 0)
        {
          main_index_dict.at(gnss_sat_id) = observed_gnss_index;
          pre_main_observing_ch = now_main_observing_ch;
          AllocateToCh(gnss_sat_id, now_main_observing_ch, main_free_ch);
        }
        double gnss_clock = gnss_satellites_.Get_true_info().GetSatelliteClock(gnss_sat_id); // これはbiasらしい
        double l1_pseudo_range = gnss_satellites_.GetPseudoRangeECI(gnss_sat_id, sat_position_i, sat_clock_true, L1_frequency);
        double l2_pseudo_range = gnss_satellites_.GetPseudoRangeECI(gnss_sat_id, sat_position_i, sat_clock_true, L2_frequency);
        auto l1_carrier_phase = gnss_satellites_.GetCarrierPhaseECI(gnss_sat_id, sat_position_i, sat_clock_true, L1_frequency);
        auto l2_carrier_phase = gnss_satellites_.GetCarrierPhaseECI(gnss_sat_id, sat_position_i, sat_clock_true, L2_frequency);

        double ionfree_range = (pow(L1_frequency/L2_frequency, 2.0)*l1_pseudo_range - l2_pseudo_range)/(pow(L1_frequency/L2_frequency, 2.0) - 1);
        double ionfree_phase = (pow(L1_frequency/L2_frequency, 2.0)*L1_lambda*(l1_carrier_phase.first + l1_carrier_phase.second) - L2_lambda*(l2_carrier_phase.first + l2_carrier_phase.second))/(pow(L1_frequency/L2_frequency, 2.0) - 1);

        gnss_true.observable_gnss_sat_id.push_back(gnss_sat_id);
        gnss_true.gnss_satellites_position.push_back(gnss_position);
        gnss_true.gnss_clock.push_back(gnss_clock);
        gnss_true.L1_pseudo_range.push_back(l1_pseudo_range);
        gnss_true.L2_pseudo_range.push_back(l2_pseudo_range);
        gnss_true.L1_carrier_phase.push_back(l1_carrier_phase);
        gnss_true.L2_carrier_phase.push_back(l2_carrier_phase);

        gnss_true.ionfree_pseudo_range.push_back(ionfree_range);
        gnss_true.ionfree_carrier_phase.push_back(ionfree_phase);

        // 観測情報の方には観測誤差を混ぜる
        std::normal_distribution<> pseudo_range_noise(0.0, pseudo_sigma);
        std::normal_distribution<> carrier_phase_noise(0.0, carrier_sigma);

        //estimateに使う方の情報
        gnss_position = gnss_satellites_.GetSatellitePositionEci(i);
        gnss_clock = gnss_satellites_.GetSatelliteClock(i);
        
        // add measurement error
        l1_pseudo_range += pseudo_range_noise(mt);
        l2_pseudo_range += pseudo_range_noise(mt);
        l1_carrier_phase.first += carrier_phase_noise(mt)/L1_lambda;
        l2_carrier_phase.first += carrier_phase_noise(mt)/L2_lambda;

        ionfree_range = (pow(L1_frequency/L2_frequency, 2.0)*l1_pseudo_range - l2_pseudo_range)/(pow(L1_frequency/L2_frequency, 2.0) - 1);
        ionfree_phase = L2_lambda*(L1_frequency/L2_frequency* (l1_carrier_phase.first + l1_carrier_phase.second) - (l2_carrier_phase.first + l2_carrier_phase.second))/(pow(L1_frequency/L2_frequency, 2.0) - 1);

        gnss_observed_values.observable_gnss_sat_id.push_back(gnss_sat_id);
        gnss_observed_values.gnss_satellites_position.push_back(gnss_position);
        gnss_observed_values.gnss_clock.push_back(gnss_clock);
        gnss_observed_values.L1_pseudo_range.push_back(l1_pseudo_range);
        gnss_observed_values.L2_pseudo_range.push_back(l2_pseudo_range);
        gnss_observed_values.L1_carrier_phase.push_back(l1_carrier_phase);
        gnss_observed_values.L2_carrier_phase.push_back(l2_carrier_phase);

        gnss_observed_values.ionfree_pseudo_range.push_back(ionfree_range);
        gnss_observed_values.ionfree_carrier_phase.push_back(ionfree_phase);
    }

    if(!gnss_observed_values.check_normal()){
        cout << "gnss_observed_values is somthing wrong" << endl;
        abort();
    }

    return;
}

// ionfreeの計算をしている．
void PBD_dgps::ProcessGnssObservation(GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, const int sat_id)
{
    int observed_gnss_index = 0;
    for(int i = 0; i < num_of_gnss_satellites; ++i)
    {
        if(observe_info.at(sat_id).pre_observed_status.at(i) == true && observe_info.at(sat_id).now_observed_status.at(i) == false)
        {
            L1_bias[sat_id].at(i) = 0.0;
            L2_bias[sat_id].at(i) = 0.0;
        }
        else if(observe_info.at(sat_id).pre_observed_status.at(i) == false && observe_info.at(sat_id).now_observed_status.at(i) == true)
        {
            // (first + second)*lambda から真の距離引いてそこからbias求める．ここから！！！！！！！！！！！！！！！！！！！！！ そのまま真の距離引いたら0になるからここでの真の距離は時刻を使う．時刻の精度以下に埋もれる部分が整数不定性として出てくる？伝搬時間も必要やん．
            L1_bias[sat_id].at(i) = gnss_true.L1_carrier_phase.at(observed_gnss_index).second; // これじゃダメ．あまり分のbiasを求めないと．あんま関係ない気がするので後で対応する．
            // これがどのchに対応しているかはわかっている．
            L2_bias[sat_id].at(i) = gnss_true.L2_carrier_phase.at(observed_gnss_index).second;
        }
        // 今はL1
        UpdateTrueBias(L1_bias, i, L1_lambda);
        if(observe_info.at(sat_id).now_observed_status.at(i)) ++observed_gnss_index;
    }

    for(int i = 0;i < gnss_true.observable_gnss_sat_id.size();++i)
    {
        int gnss_sat_id = gnss_true.observable_gnss_sat_id.at(i);

        auto L1_observed = gnss_true.L1_carrier_phase.at(i);
        L1_observed.first += L1_observed.second - L1_bias[sat_id].at(gnss_sat_id); // 何がしたいのか？

        auto L2_observed = gnss_true.L2_carrier_phase.at(i);
        L2_observed.first += L2_observed.second - L2_bias[sat_id].at(gnss_sat_id);

        double ionfree_phase = (pow(L1_frequency/L2_frequency, 2.0)*L1_lambda*L1_observed.first - L2_lambda*L2_observed.first)/(pow(L1_frequency/L2_frequency, 2.0) - 1);
        gnss_true.ionfree_carrier_phase.push_back(ionfree_phase);
    }

    for(int i = 0;i < gnss_observed_values.observable_gnss_sat_id.size();++i)
    {
        int gnss_sat_id = gnss_observed_values.observable_gnss_sat_id.at(i);

        auto L1_observed = gnss_observed_values.L1_carrier_phase.at(i);
        L1_observed.first += L1_observed.second - L1_bias[sat_id].at(gnss_sat_id);

        auto L2_observed = gnss_observed_values.L2_carrier_phase.at(i);
        L2_observed.first += L2_observed.second - L2_bias[sat_id].at(gnss_sat_id);

        double ionfree_phase = (pow(L1_frequency/L2_frequency, 2.0)*L1_lambda*L1_observed.first - L2_lambda*L2_observed.first)/(pow(L1_frequency/L2_frequency, 2.0) - 1);
        gnss_observed_values.ionfree_carrier_phase.push_back(ionfree_phase);
    }
    return;
}

// ここがバグっている？
void PBD_dgps::UpdateTrueBias(vector<vector<double>> bias, const int gnss_sat_id, const double lambda)
{
  if (now_main_observing_ch.count(gnss_sat_id))
  {
    int ch = now_main_observing_ch[gnss_sat_id];
    true_bias_main(ch) = bias[0].at(gnss_sat_id)*lambda;
  }
  if (now_common_observing_ch.count(gnss_sat_id))
  {
    int ch = now_common_observing_ch[gnss_sat_id];
    true_bias_target(ch) = (bias[1].at(gnss_sat_id) - bias[0].at(gnss_sat_id))*lambda;
  }
};

// sat_idはLEO衛星のこと
void PBD_dgps::FindCommonObservedGnss(const std::pair<int, int> sat_id_pair)
{
  // 初期化
  common_observed_gnss_sat_id.clear();
  common_observed_status.assign(num_of_gnss_satellites, false);
  const int main_sat_id = sat_id_pair.first;
  const int target_sat_id = sat_id_pair.second;
  int common_index = 0;
  for (int i = 0; i < observe_info.at(main_sat_id).now_observed_gnss_sat_id.size(); ++i)
  {
    for (int j = 0; j < observe_info.at(target_sat_id).now_observed_gnss_sat_id.size(); ++j)
    {
      int gnss_sat_id = observe_info.at(main_sat_id).now_observed_gnss_sat_id.at(i);
      // どっかで複数衛星にも拡張
      if (gnss_sat_id == observe_info.at(target_sat_id).now_observed_gnss_sat_id.at(j))
      {
        common_observed_gnss_sat_id.push_back(gnss_sat_id);
        common_observed_status.at(gnss_sat_id) = true;
        common_index_dict.at(gnss_sat_id) = common_index; // このindexはほんまに必要なのか？
        ++common_index;
        pre_common_observing_ch = now_common_observing_ch;
        AllocateToCh(gnss_sat_id, now_common_observing_ch, common_free_ch);
        break;
      }
      pre_common_observing_ch = now_common_observing_ch;
      RemoveFromCh(gnss_sat_id, now_common_observing_ch, common_free_ch);
    }
  }
}

void PBD_dgps::AllocateToCh(const int gnss_sat_id, std::map<const int, int>& observing_ch, vector<int>& free_ch)
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

void PBD_dgps::RemoveFromCh(const int gnss_sat_id, std::map<const int, int>& observing_ch, vector<int>& free_ch)
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

// mainとtargetに関して分けてするのがいいのでは？
void PBD_dgps::UpdateBiasForm(const int sat_id, EstimatedVariables& x_est)// LEO衛星の数が増えたときは衛星ごとにこのクラスのインスタンスが生成される？ので一旦これで行く
{
  //観測する衛星同じだったら飛ばしていい
  if (CheckVectorEqual(observe_info.at(sat_id).pre_observed_gnss_sat_id, observe_info.at(sat_id).now_observed_gnss_sat_id))
  {
    for (int i = 0; i < num_of_gnss_satellites; ++i)
    {
      observe_info.at(sat_id).now_observed_status.at(i) = false;
    }
    return;
  }
  int n = observe_info.at(sat_id). now_observed_gnss_sat_id.size();
  int n_pre = observe_info.at(sat_id).pre_observed_gnss_sat_id.size();

  // 多分ここは使っていない
  // vector<int> pre_index_from_sat_id(num_of_gnss_satellites, -1);
  // int index = 0;
  // for (auto gnss_id : pre_observed_gnss_sat_id.at(sat_id)) {
  //   pre_index_from_sat_id.at(gnss_id) = index;
  //   ++index;
  // }

  // vector<int> now_index_from_sat_it(num_of_gnss_satellites, -1);
  // index = 0;
  // for (auto gnss_id : now_observed_gnss_sat_id.at(sat_id))
  // {
  //   now_index_from_sat_it.at(gnss_id) = index;
  //   ++index;
  // }

  int pre_index = 0; //index is not id
  int now_index = 0; //index is not id

  const Eigen::VectorXd& pre_estimated_bias = x_est.bias;
  Eigen::MatrixXd pre_M = M;

  /*
  double geo_ure_sigma = sqrt(pow(L1_frequency / L2_frequency, 4.0) + 1.0) * pseudo_sigma / (pow(L1_frequency / L2_frequency, 2.0) - 1.0); //ionfree_pseudo_sigma
  */
  // i = gnss_sat_id
  for (int i = 0; i < num_of_gnss_satellites; ++i)
  {
    if (observe_info.at(sat_id).pre_observed_status.at(i) == false && observe_info.at(sat_id).now_observed_status.at(i) == false) continue;
    // 見えなくなったとき
    else if (observe_info.at(sat_id).pre_observed_status.at(i) == true && observe_info.at(sat_id).now_observed_status.at(i) == false)
    {
      // if (pre_index != pre_index_from_sat_id.at(i)) {
      //   cout << "pre index something is wrong" << endl;
      //   abort();
      // }
      x_est.bias(pre_index) = 0.0;
      // あと全部初期化するんではなくて，Mの成分がどこが変化しうるのかを考えて計算した方がいい．
      M.block(0, num_of_single_status + pre_index, state_dimension, 1) = Eigen::MatrixXd::Zero(state_dimension, 1);
      M.block(num_of_single_status + pre_index, 0, state_dimension, 1) = Eigen::MatrixXd::Zero(1, state_dimension);
      M(num_of_single_status + pre_index, num_of_single_status + pre_index) = 0.0; //pow(sigma_N_process, 2.0); // 対角だけ残す．1とかでいいのか？
      ++pre_index;
    }
    else if (observe_info.at(sat_id).pre_observed_status.at(i) == false && observe_info.at(sat_id).now_observed_status.at(i) == true)
    {
      // if (now_index != now_index_from_sat_it.at(i))
      // {
      //   cout << "now index something is wrong 1" << endl;
      //   abort();
      // }
      std::normal_distribution<> N_dist(0.0, sigma_N_process);
      x_est.bias(now_index) = N_dist(mt); // こいつもGauss Makov過程にした方がいいらしい?
      M(num_of_single_status + now_index, num_of_single_status + now_index) = pow(sigma_N_process, 2.0); // 列と行に関して0にしなくていい?
      ++now_index;
    }
    else if (observe_info.at(sat_id).pre_observed_status.at(i) == true && observe_info.at(sat_id).now_observed_status.at(i) == true)
    {
      x_est.bias(now_index) = pre_estimated_bias(pre_index);
      // あと全部初期化するんではなくて，Mの成分がどこが変化しうるのかを考えて計算した方がいい．
      M.block(0, num_of_single_status + now_index, state_dimension, 1) = pre_M.block(0, num_of_single_status + pre_index, state_dimension, 1);
      M.block(num_of_single_status + now_index, 0, state_dimension, 1) = pre_M.block(num_of_single_status + pre_index, 0, state_dimension, 1);
      ++pre_index; ++now_index;
      /*
      //jは今存在するところしか走らない
      for (int j = 0; j < now_observed_gnss_sat_id.at(sat_id).size(); ++j)
      {
        if (j == now_index) break; //自分に追いついたのでbreak
        int j_now_gnss_sat_id = now_observed_gnss_sat_id.at(sat_id).at(j); //向こうの今のsat_id
        // いまいち何してんのかわからん<-MのNに関する成分をコピっている．関数化したい．
        if (pre_index_from_sat_id.at(j_now_gnss_sat_id) == -1) continue;
        int j_pre_index = pre_index_from_sat_id.at(j_now_gnss_sat_id);

        // 対角じゃないところ．
        M(j + num_of_single_status + M_offset, now_index + num_of_single_status + M_offset) = pre_M(j_pre_index + num_of_single_status + pre_M_offset, pre_index + num_of_single_status + pre_M_offset);
        M(now_index + num_of_single_status + M_offset, j + num_of_single_status + M_offset) = pre_M(pre_index + num_of_single_status + pre_M_offset, j_pre_index + num_of_single_status + pre_M_offset);
      }

      //対角成分引継ぎ
      M(now_index + num_of_single_status + M_offset, now_index + num_of_single_status + M_offset) = pre_M(pre_index + num_of_single_status + pre_M_offset, pre_index + num_of_single_status + pre_M_offset);

      if (pre_index != pre_index_from_sat_id.at(i)) {
        cout << "pre index something is wrong" << endl;
        abort();
      }
      if (now_index != now_index_from_sat_it.at(i)) {
        cout << "now index something is wrong 2" << endl;
        cout << i << endl;
        for (int j = 0; j < now_observed_gnss_sat_id.at(sat_id).size(); ++j) {
          cout << j << " " << now_observed_gnss_sat_id.at(sat_id).at(j) << endl;
        }
        cout << now_index << endl;
        abort();
      }
      */
    }
  }
  /*
  if (abs(M.determinant()) < 10e-13)
  {
    cout << "M matirx is singular" << endl;
    abort();
  }
  */
  //cout << M << endl;
}

void PBD_dgps::SetBiasToObservation(const int sat_id, EstimatedVariables& x_est, GnssObservedValues gnss_observed, GnssObservedValues gnss_true)
{
    // 共通衛星見つける
    FindCommonObservedGnss(std::make_pair(0, 1)); // ここではない気がする．

    Eigen::MatrixXd pre_M = M;
    int n = observe_info.at(sat_id).now_observed_gnss_sat_id.size();
    int n_pre = observe_info.at(sat_id).pre_observed_gnss_sat_id.size();
    
    // M = Eigen::MatrixXd::Zero(state_dimension, state_dimension); //0クリア <- これ必要かな？
    // bias以外のstatusはそのまま
    /*
    M.topLeftCorner(num_of_single_status, num_of_single_status) = pre_M.topLeftCorner(num_of_single_status, num_of_single_status);
    M.block(single_dimension, single_dimension, num_of_single_status, num_of_single_status) = pre_M.block(single_dimension, single_dimension, num_of_single_status, num_of_single_status);
    */

    UpdateBiasForm(sat_id, x_est);
    // update observation state info
    for (int i = 0; i < num_of_gnss_satellites; ++i)
    {
      observe_info.at(sat_id).pre_observed_status.at(i) = observe_info.at(sat_id).now_observed_status.at(i);
      observe_info.at(sat_id).now_observed_status.at(i) = false;
    }
    observe_info.at(sat_id).pre_observed_gnss_sat_id.clear();
    for (auto x : observe_info.at(sat_id).now_observed_gnss_sat_id) observe_info.at(sat_id).pre_observed_gnss_sat_id.push_back(x);
    return;
}

bool PBD_dgps::CheckCanSeeSatellite(const libra::Vector<3> satellite_position, const libra::Vector<3> gnss_position) const
{
    double angle_rad = angle(satellite_position, gnss_position - satellite_position);
    if(angle_rad < M_PI/2.0 - mask_angle) return true;
    else return false;
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


template <typename T> bool PBD_dgps::CheckVectorEqual(const vector<T>& a, const vector<T>& b)
{
    if(a.size() != b.size()) return false;
    for(int i = 0;i < a.size();++i){
        if(a.at(i) != b.at(i)) return false;
    }

    return true;
}