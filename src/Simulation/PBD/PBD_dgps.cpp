#include <iomanip>
#include "PBD_dgps.h"
#include "PBD_const.h"

// outputを変えるときは"result.csv"を変更する．せめてパスは変えたい．
PBD_dgps::PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Orbit& main_orbit, const Orbit& target_orbit):mt(42), sat_main_clock_true(0.0), sat_target_clock_true(0.0), step_time(sim_time_.GetStepSec()), ofs("result_new.csv"), num_of_gnss_sastellites(gnss_satellites_.GetNumOfSatellites()), main_orbit_(main_orbit), target_orbit_(target_orbit)
{
    //初期化
    estimated_status = Eigen::VectorXd::Zero(4);
    libra::Vector<3> position_main = main_orbit_.GetSatPosition_i();
    for(int i = 0;i < 3;++i) estimated_status(i) = position_main[i];
    libra::Vector<3> velocity_main = main_orbit_.GetSatVelocity_i();
    for(int i = 0;i < 3;++i) estimated_velocity(i) = velocity_main[i];
    estimated_acc = Eigen::VectorXd::Zero(3);

    estimated_differential_status = Eigen::VectorXd::Zero(4);
    libra::Vector<3> position_target = target_orbit_.GetSatPosition_i();
    for (int i = 0; i < 3; ++i) estimated_differential_status(i) = position_target[i] - position_main[i];
    libra::Vector<3> velocity_target = target_orbit_.GetSatVelocity_i();
    for (int i = 0; i < 3; ++i) estimated_differential_velocity(i) = velocity_target[i] - velocity_main[i];
    estimated_differential_acc = Eigen::VectorXd::Zero(3);

    // 初期分散
    std::normal_distribution<> position_dist(0.0,sigma_r_ini);
    std::normal_distribution<> clock_dist(0.0, sigma_cdt_ini);
    std::normal_distribution<> velocity_dist(0.0, sigma_v_ini);

    // Vって名前は変えた方がいいかも
    Eigen::VectorXd V = Eigen::VectorXd::Constant(num_of_status, 0);

    // A-priori 
    for(int i = 0;i < 3;++i) V(i) = pow(sigma_r_ini, 2.0); // position
    V(3) = pow(sigma_cdt_ini, 2.0); // clock
    for(int i = 4;i < 7;++i) V(i) = pow(sigma_v_ini, 2.0); // velocity
    for (int i = 7; i < 10; ++i) V(i) = pow(sigma_acc_r_ini, 2.0); // acceleration RTNになっていないので後ほど修正　
    for(int i = 7;i < 10;++i) V(i) = pow(sigma_acc_r_ini, 2.0);
    // Nについても入れる．
    // 以下がdifferentialの状態量
    for (int i = num_of_single_status; i < num_of_single_status + 3; ++i) V(i) = pow(sigma_r_ini, 2.0);
    V(num_of_single_status + 3) = pow(sigma_cdt_ini, 2.0);
    for (int i = num_of_single_status + 4; i < num_of_single_status + 7; ++i) V(i) = pow(sigma_v_ini, 2.0);
    for (int i = num_of_single_status + 7; i < num_of_single_status + 10; ++i) V(i) = pow(sigma_acc_r_ini, 2.0); // ここもRTNになってないので後ほど修正

    M = V.asDiagonal(); // 誤差共分散行列M

    // 初期位置はガウシアンからサンプリング．mtは多分，乱数のシード
    for(int i = 0; i < 3; ++i) estimated_status(i) += position_dist(mt);
    for(int i = 0; i < 3; ++i) estimated_velocity(i) += velocity_dist(mt);
    for(int i = 0; i < 3; ++i) estimated_differential_status(i) += position_dist(mt);
    for(int i = 0; i < 3; ++i) estimated_differential_velocity(i) += velocity_dist(mt);

    for (int i = 0; i < 2; ++i) pre_observed_vector.at(i).assign(num_of_gnss_sastellites, false);
    for (int i = 0; i < 2; ++i) now_observed_vector.at(i).assign(num_of_gnss_sastellites, false);
    common_observed_vector.at(0).assign(num_of_gnss_sastellites, false);
    for (int i = 0; i < 2; ++i) first_L1_bias.at(i).assign(num_of_gnss_sastellites, 0.0);
    for (int i = 0; i < 2; ++i) first_L2_bias.at(i).assign(num_of_gnss_sastellites, 0.0);

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
    sat_main_clock_true = clock_dist(mt); // main
    sat_target_clock_true = clock_dist(mt); // target seadsも変えた方がいいかも．というかdefault_random_engineを使うのがよさそう．

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
        ProcessGnssObservation(gnss_observed_values_main, gnss_true_main, 0);
        GnssObservedValues gnss_observed_values_target;
        GnssObservedValues gnss_true_target;
        GetGnssPositionObservation(gnss_satellites_, target_orbit_, 1, gnss_observed_values_target, gnss_true_target, sat_target_clock_true);
        ProcessGnssObservation(gnss_observed_values_target, gnss_true_target, 1);
        std::pair<GnssObservedValues, GnssObservedValues> gnss_observed_pair = std::make_pair(gnss_observed_values_main, gnss_observed_values_target);
        std::pair<GnssObservedValues, GnssObservedValues> gnss_true_pair = std::make_pair(gnss_true_main, gnss_true_target);
        resize_Matrix(gnss_observed_pair, gnss_true_pair);
        KalmanFilter(gnss_observed_values_main, gnss_observed_values_target);
    }

    //log output
    if(abs(elapsed_time - tmp_log*log_step_time) < step_time/2.0){
        libra::Vector<3> sat_position = main_orbit_.GetSatPosition_i();
        libra::Vector<3> sat_velocity = main_orbit_.GetSatVelocity_i();
        for(int i = 0;i < 3;++i) ofs << fixed << setprecision(30) << sat_position[i] << ",";
        ofs << fixed << setprecision(30) << sat_main_clock_true << ",";
        for(int i = 0;i < 3;++i) ofs << fixed << setprecision(30) << sat_velocity[i] << ",";
        libra::Vector<3> sat_position_target = target_orbit_.GetSatPosition_i();
        libra::Vector<3> sat_velocity_target = target_orbit_.GetSatVelocity_i();
        for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << sat_position_target[i] << ",";
        ofs << fixed << setprecision(30) << sat_target_clock_true << ","; // ?
        for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << sat_velocity_target[i] << ",";
        for(int i = 0;i < 4;++i) ofs << fixed << setprecision(30) << estimated_status(i) << ",";
        for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << estimated_velocity(i) << ",";
        for(int i = 0;i < 3;++i) ofs << fixed << setprecision(30) << estimated_acc(i) << ",";
        for (int i = 0; i < 4; ++i) ofs << fixed << setprecision(30) << estimated_differential_status(i) << ",";
        for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << estimated_differential_velocity(i) << ",";
        for (int i = 0; i < 3; ++i) ofs << fixed << setprecision(30) << estimated_differential_acc(i) << ",";
        for (int i = 0; i < num_of_single_status; ++i) ofs << fixed << setprecision(30) << M(i, i) << ",";
        int main_observed_gnss_num = now_observed_gnss_sat_id.at(0).size();
        for(int i = num_of_single_status + main_observed_gnss_num;i < num_of_status + main_observed_gnss_num; ++i) ofs << fixed << setprecision(30) << M(i, i) << ",";
        // ここに整数不定性も併せてログに入れたいが，そもそも観測衛星数が変化するので全英成分をまとめてログに入れるようにしないといけない．
        // ここは何を残しているのか確認．
        int ans = 0;
        for(int i = 0;i < num_of_gnss_sastellites;++i){
            auto gnss_position = gnss_satellites_.GetSatellitePositionEci(i);
            if(CheckCanSeeSatellite(sat_position, gnss_position)) ++ans;
        }
        ofs << ans << endl;
    }
    
    return;
}

void PBD_dgps::OrbitPropagation()
{
    //RK4
    Eigen::Vector3d position_main = estimated_status.topRows(3);
    Eigen::Vector3d velocity_main = estimated_velocity;
    Eigen::Vector3d acceleration_main = estimated_acc;
    Eigen::Vector3d position_target = estimated_differential_status.topRows(3) + position_main;
    Eigen::Vector3d velocity_target = estimated_differential_velocity + velocity_main;
    Eigen::Vector3d acceleration_target = estimated_differential_acc + acceleration_main;

    // ここはcoreの機能を使うように修正
    vector<Eigen::Vector3d> pos_vel_main = RK4(position_main, velocity_main, acceleration_main);
    position_main = pos_vel_main[0];
    velocity_main = pos_vel_main[1];

    vector<Eigen::Vector3d> pos_vel_target = RK4(position_target, velocity_target, acceleration_target);
    position_target = pos_vel_target[0];
    velocity_target = pos_vel_target[1];

    for(int i = 0;i < 3;++i) estimated_status(i) = position_main(i);
    estimated_velocity = velocity_main;
    for (int i = 0; i < 3; ++i) estimated_differential_status(i) = position_target(i) - position_main(i);
    estimated_differential_velocity = velocity_target - velocity_main;

    M = update_M_matrix(position_main, velocity_main, acceleration_main, position_target - position_main, estimated_differential_velocity, estimated_differential_acc);
}

// Orbitの更新を使えるように修正
vector<Eigen::Vector3d> PBD_dgps::RK4(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration)
{
  Eigen::Vector3d k0 = position_differential(velocity);
  Eigen::Vector3d l0 = velocity_differential(position, velocity, acceleration);

  Eigen::Vector3d tmp_position = position + k0 * step_time / 2.0;
  Eigen::Vector3d tmp_velocity = velocity + l0 * step_time / 2.0;
  Eigen::Vector3d k1 = position_differential(tmp_velocity);
  Eigen::Vector3d l1 = velocity_differential(tmp_position, tmp_velocity, acceleration);

  tmp_position = position + k1 * step_time / 2.0;
  tmp_velocity = velocity + l1 * step_time / 2.0;
  Eigen::Vector3d k2 = position_differential(tmp_velocity);
  Eigen::Vector3d l2 = velocity_differential(tmp_position, tmp_velocity, acceleration);

  tmp_position = position + k2 * step_time;
  tmp_velocity = velocity + l2 * step_time;
  Eigen::Vector3d k3 = position_differential(tmp_velocity);
  Eigen::Vector3d l3 = velocity_differential(tmp_position, tmp_velocity, acceleration);

  vector<Eigen::Vector3d> position_velocoty = 
  { position + step_time * (k0 + 2.0 * k1 + 2.0 * k2 + k3) / 6.0, 
    velocity + step_time * (l0 + 2.0 * l1 + 2.0 * l2 + l3) / 6.0 };
  return position_velocoty;
}

Eigen::Vector3d PBD_dgps::position_differential(const Eigen::Vector3d& velocity) const
{
    return velocity;
}

// 単位が違うのでそこの変換が必要　
Eigen::Vector3d PBD_dgps::velocity_differential(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration) const
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

    return all_acceleration;
}


Eigen::MatrixXd PBD_dgps::update_M_matrix(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration, const Eigen::Vector3d& position_difference, const Eigen::Vector3d& velocity_difference, const Eigen::Vector3d& acceleration_difference)
{
  int n_main = estimated_bias.size(); //mainの見るGNSS数
  int n_common = estimated_differential_bias.size(); //共通のGNSS数

  // Dinamics and Kinematics model
  Eigen::MatrixXd A_jacobi = calculate_A_matrix(position, velocity, acceleration, position_difference, velocity_difference, acceleration_difference); //Jacobi行列
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_of_status + n_main + n_common, num_of_status + n_main + n_common); 
  A.topLeftCorner(num_of_single_status, num_of_single_status) = A_jacobi.topLeftCorner(num_of_single_status, num_of_single_status);
  A.block(num_of_single_status + n_main, num_of_single_status + n_main, num_of_single_status, num_of_single_status) = A_jacobi.bottomRightCorner(num_of_single_status, num_of_single_status); 
  // STM
  Eigen::MatrixXd Phi = Eigen::MatrixXd::Identity(num_of_status + n_main + n_common, num_of_status + n_main + n_common);
  // clock
  Phi(3, 3) = 1.0;
  Phi(num_of_single_status + n_main + 3, num_of_single_status + n_main + 3) = 1.0;
  Phi += step_time * A;
  // AからPhiを求めるというこの構成も変えた方がいいのかもしれない．経験加速度に関してはdtを掛けない気がするし．
  // a
  Phi.block(7, 7, 3, 3) = calculate_Phi_a(step_time);
  Phi.block(num_of_single_status + n_main + 7, num_of_single_status + n_main + 7, 3, 3) = calculate_Phi_a(step_time);

  // ここもaについて修正．
  Eigen::MatrixXd Q = calculate_Q_matrix(n_main, n_common);
  
  // Gamma = BQB^t
  Eigen::MatrixXd Gamma = pow(step_time, 2.0) * Q;
  // Gamma(3, 3) = pow(clock_sigma, 2.0); // ここも推定量ならこれはおかしい？
  // Gamma(num_of_single_status + n_main + 3, num_of_single_status + n_main + 3) = pow(clock_sigma, 2.0);

  Eigen::MatrixXd res = Phi * M * Phi.transpose() + Gamma;

  return res;
}

// Jacobi or STM for pod
Eigen::MatrixXd PBD_dgps::calculate_A_matrix(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration) const
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

Eigen::MatrixXd PBD_dgps::calculate_A_matrix(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration, const Eigen::Vector3d& position_difference, const Eigen::Vector3d& velocity_difference, const Eigen::Vector3d& acceleration_difference) const
{
  double r = position.norm(); // [m]
  double v = velocity.norm(); // [m/s]
  double a = acceleration.norm(); // [nm/s]

  double x = position(0); double y = position(1); double z = position(2);
  double vx = velocity(0); double vy = velocity(1); double vz = velocity(2);

  double J2_coefficient = 3.0 / 2.0 * mu_const * J2_const * Earth_Radius * Earth_Radius;

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_of_status, num_of_status);
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
  
  // Kinematics model
  // (dr, dv)
  A(num_of_single_status, num_of_single_status + 4)     = 1.0;
  A(num_of_single_status + 1, num_of_single_status + 5) = 1.0;
  A(num_of_single_status + 2, num_of_single_status + 6) = 1.0;
  // (dv, da)
  A(num_of_single_status + 4, num_of_single_status + 7) = 1.0 * 1e-6;
  A(num_of_single_status + 5, num_of_single_status + 8) = 1.0 * 1e-6;
  A(num_of_single_status + 6, num_of_single_status + 9) = 1.0 * 1e-6;

  return A;
}

// Process noiseによる共分散行列へのノイズ分
Eigen::MatrixXd PBD_dgps::calculate_Q_matrix(const int n_main, const int n_common)
{
  const int n = n_main + n_common;
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(n + num_of_status, n + num_of_status);
  for(int i = 0;i < 3;++i) Q(i, i) = pow(sigma_r_process, 2.0);
  Q(3, 3) = pow(sigma_cdt_process, 2.0);
  for(int i = 4;i < 7;++i) Q(i, i) = pow(sigma_v_process, 2.0);
  for(int i = 7;i < 10;++i) Q(i, i) = pow(sigma_acc_r_process, 2.0);
  for (int i = num_of_single_status; i < n_main + num_of_single_status; ++i) Q(i, i) = pow(sigma_N_process, 2.0);
  int single_offset = num_of_single_status + n_main;
  for (int i = num_of_single_status; i < num_of_single_status + 3; ++i) Q(i, i) = pow(sigma_r_process, 2.0);
  Q(single_offset + 3, single_offset + 3) = pow(sigma_cdt_process, 2.0);
  for (int i = num_of_single_status + 4; i < num_of_single_status + 7; ++i) Q(i, i) = pow(sigma_v_process, 2.0);
  for(int i = single_offset + 7; i < single_offset + num_of_single_status; ++i) Q(i, i) = pow(sigma_acc_r_process, 2.0);
  for (int i = num_of_single_status + single_offset; i < num_of_status + n_main + n_common; ++i) Q(i, i) = pow(sigma_N_process, 2.0);

  return Q;
}

Eigen::MatrixXd PBD_dgps::calculate_Phi_a(const double dt)
{
  // ここが各軸に応じて変わる部分なのか．
  const double tau = 900; // [s] <- constに移すべきか？
  double phi = exp(-dt / tau);
  Eigen::MatrixXd Phi_a = Eigen::MatrixXd::Identity(3, 3);
  Phi_a *= phi;
  return Phi_a;
};

void PBD_dgps::KalmanFilter(const GnssObservedValues& gnss_observed_main, const GnssObservedValues& gnss_observed_target)
{
    // gnss_observed_hoge の観測される時間はどうなっているのか？
    int n_main = gnss_observed_main.can_see_satellites_sat_id.size();
    int n_common = common_observed_gnss_sat_id[0].size();
    if(n_main != estimated_bias.size()){
        cout << "estimated bias and observe is something wrong" << endl;
        abort();
    }
    if (n_common != estimated_differential_bias.size()) {
      cout << "estimated differential bias and observe is something wrong" << endl;
      abort();
    }
    
    Eigen::VectorXd z = Eigen::VectorXd(2*n_main + n_common); //観測ベクトル
    Eigen::VectorXd h_x = Eigen::VectorXd(2*n_main + n_common); // 観測モデル行列
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(2*n_main + n_common, num_of_status + n_main + n_common); //観測行列（dhx/dx）

    for(int i = 0;i < n_main;++i){
        auto gnss_position = gnss_observed_main.gnss_satellites_position.at(i);
        double gnss_clock = gnss_observed_main.gnss_clock.at(i);
        // ionfreeでいいのか？←今後電離圏モデルの追加も必要になる．
        double pseudo_range = gnss_observed_main.ionfree_pseudo_range.at(i);
        double carrier_phase = gnss_observed_main.ionfree_carrier_phase.at(i);
        // 共通衛星を探す．
        unsigned int gnss_sat_id = gnss_observed_main.can_see_satellites_sat_id.at(i);
        double gnss_clock_target = 0.0;
        double carrier_phase_target = 0.0;
        for (int j = 0; j < gnss_observed_target.can_see_satellites_sat_id.size(); ++j)
        {
          // 共通衛星を発見
          if (gnss_sat_id == gnss_observed_target.can_see_satellites_sat_id.at(j))
          {
            gnss_clock_target = gnss_observed_target.gnss_clock.at(j);
            carrier_phase_target = gnss_observed_target.ionfree_carrier_phase.at(j);
            break;
          }
        }
        // 仮にここで生成する．
        std::normal_distribution<> pseudo_range_noise(0.0, pseudo_sigma);
        std::normal_distribution<> carrier_phase_noise(0.0, carrier_sigma);

        double geometric_range = calculate_geometric_range(estimated_status, gnss_position); // TODO: ここのGNSSの位置情報とかも実際には真値とずれた推定値を使わないといけない！
        double pseudo_range_model = calculate_pseudo_range(estimated_status, gnss_position, gnss_clock) + pseudo_range_noise(mt);
        double ionfree_lambda = L1_lambda; // 要修正
        double phase_range_model = calculate_carrier_phase(estimated_status, gnss_position, gnss_clock, estimated_bias(i), ionfree_lambda) + carrier_phase_noise(mt); // pseudo_range_model + estimated_bias(i);
        double pseudo_range_target_model = calculate_pseudo_range(estimated_differential_status + estimated_status, gnss_position, gnss_clock_target) + pseudo_range_noise(mt);
        double phase_range_target_model = calculate_carrier_phase(estimated_differential_status + estimated_status, gnss_position, gnss_clock, estimated_differential_bias(common_index_dict.at(gnss_sat_id)) + estimated_bias(i), ionfree_lambda) + carrier_phase_noise(mt); // pseudo_range_target_calc + estimated_differential_bias(common_index_dict.at(gnss_sat_id)) + estimated_bias(i);

        // 真の観測量
        z(i) = pseudo_range;
        z(i + n_main) = carrier_phase;
        // hは観測モデルに基づいて推定した観測量のベクトル
        h_x(i) = pseudo_range_model;
        h_x(i + n_main) = phase_range_model;
        static int k = 0; // common 衛星のID

        for(int j = 0; j < 3; ++j){
          // forじゃなくて行列計算したい．
            H(i, j) = (estimated_status(j) - gnss_position[j])/geometric_range;
            H(i + n_main, j) = (estimated_status(j) - gnss_position[j])/geometric_range;
            // SDCP
            if (gnss_clock_target != 0.0 && carrier_phase_target != 0.0) // これじゃなくて共通衛星見つけたフラグ的なものを作ったほうがいい？
            {
              double sdcp = carrier_phase_target - carrier_phase; // ホンマはここでSDCPに入りうる観測ノイズを求めて入れるべき．(sqrt(2)*epsilon_cp)
              z(2 * n_main + k) = sdcp;
              h_x(2 * n_main + k) = phase_range_target_model - phase_range_model;
              // SDCPの場合は視線方向ベクトルでいい
              H(2 * n_main + k, j + num_of_single_status + n_main) = (estimated_status(j) - gnss_position[j]) / geometric_range;
            }
        }
        k++;
        if (k >= n_common - 1 || i == n_main-1) k = 0;
        // clock
        H(i, 3) = 1.0;
        H(i + n_main, 3) = 1.0;

        // N
        H(i + n_main, i + num_of_single_status) = 1.0;
    }
    
    for (int i = 0; i < n_common; ++i)
    {
      //differential clock
      H(2 * n_main + i, num_of_single_status + n_main + 3) = 1.0;
      //dN
      H(2 * n_main + i, i + num_of_status + n_main) = 1.0;
    }

    double ionfree_pseudo_sigma = sqrt(pow(L1_frequency/L2_frequency, 4.0) + 1.0)*pseudo_sigma/(pow(L1_frequency/L2_frequency, 2.0) - 1.0);
    double ionfree_carrier_sigma = sqrt(pow(L1_frequency/L2_frequency, 4.0) + 1.0)*carrier_sigma/(pow(L1_frequency/L2_frequency, 2.0) - 1.0);
    
    Eigen::VectorXd R_V = Eigen::VectorXd::Constant(2*n_main + n_common, pow(ionfree_pseudo_sigma, 2.0));
    Eigen::VectorXd R_V_phase = Eigen::VectorXd::Constant(n_main, pow(ionfree_carrier_sigma, 2.0));
    Eigen::VectorXd R_V_phase_difference = Eigen::VectorXd::Constant(n_common, pow(ionfree_carrier_sigma, 2.0));
    R_V.block(n_main, 0, n_main, 1) = R_V_phase;
    R_V.bottomRows(n_common) = R_V_phase_difference;
    Eigen::MatrixXd R = R_V.asDiagonal();

    //カルマンゲイン
    Eigen::MatrixXd tmp = R + H*M*H.transpose();
    Eigen::MatrixXd K = M * H.transpose() * tmp.inverse();

    Eigen::VectorXd x_predict = Eigen::VectorXd(num_of_status + n_main + n_common);
    x_predict.topRows(4) = estimated_status;
    x_predict.block(4, 0, 3, 1)  = estimated_velocity;
    x_predict.block(7, 0, 3, 1) = estimated_acc;
    x_predict.block(10, 0, n_main, 1) = estimated_bias;
    //x_predict(10) = Cd;
    x_predict.block(num_of_single_status + n_main, 0, 4, 1) = estimated_differential_status;
    x_predict.block(num_of_single_status + n_main + 4, 0, 3, 1) = estimated_differential_velocity;
    x_predict.block(num_of_single_status + n_main + 7, 0, 3, 1) = estimated_differential_acc;
    x_predict.bottomRows(n_common) = estimated_differential_bias;
    
    //h_x = H * x_predict;
    // ここの更新を一緒にしているのが間違いなのでは？update後の量が差分になっていない
    Eigen::VectorXd x_update = x_predict + K*(z - h_x);

    //更新
    estimated_status = x_update.topRows(4);
    estimated_velocity = x_update.block(4, 0, 3, 1);
    estimated_acc = x_update.block(7, 0, 3, 1);
    estimated_bias = x_update.block(10, 0, n_main, 1);
    //Cd = x_update(10);
    estimated_differential_status = x_update.block(num_of_single_status + n_main, 0, 4, 1);
    estimated_differential_velocity = x_update.block(num_of_single_status + n_main + 4, 0, 3, 1);
    estimated_differential_acc = x_update.block(num_of_single_status + n_main + 7, 0, 3, 1);
    estimated_differential_bias = x_update.bottomRows(n_common);

    //Eigen::MatrixXd tmp2 = M * H.transpose() * tmp.inverse() * H * M;
    M = M - K * H * M; // ここの更新が怪しい．
    //Eigen::MatrixXd H_trans = H.transpose();
    //M = M - H.inverse()*R*H_trans.inverse();
    
    if (!std::isfinite(estimated_status(0)))
    {
      cout << "inf or nan" << endl;
      abort();
    }
    return;
} // Mの差分成分が全く更新されず一定なのが変．

void PBD_dgps::GetGnssPositionObservation(const GnssSatellites& gnss_satellites_, const Orbit& orbit_, const int sat_id, GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, double sat_clock_true)
{
    //推定値の計算
    int num_of_gnss_satellites = gnss_satellites_.GetNumOfSatellites();
    libra::Vector<3> sat_position_i = orbit_.GetSatPosition_i();

    now_observed_gnss_sat_id[sat_id].clear(); //クリア

    for(int i = 0; i < num_of_gnss_satellites; ++i){
        //if(i == 7 || i == 23 || i == 31) continue;
        if(!gnss_satellites_.GetWhetherValid(i)) continue;
        libra::Vector<3> gnss_position = gnss_satellites_.Get_true_info().GetSatellitePositionEci(i);
        bool see_flag = CheckCanSeeSatellite(sat_position_i, gnss_position);
        // init
        main_index_dict.insert(std::make_pair(i, -1));
        common_index_dict.insert(std::make_pair(i, -1));

        if(!see_flag) continue;

        int gnss_sat_id = i;
        now_observed_vector[sat_id].at(gnss_sat_id) = true;
        now_observed_gnss_sat_id[sat_id].push_back(gnss_sat_id);
        int observed_gnss_index = now_observed_gnss_sat_id[sat_id].size() - 1;
        if (sat_id == 0) main_index_dict.at(gnss_sat_id) = observed_gnss_index;
        double gnss_clock = gnss_satellites_.Get_true_info().GetSatelliteClock(gnss_sat_id);
        double l1_pseudo_range = gnss_satellites_.GetPseudoRangeECI(gnss_sat_id, sat_position_i, sat_clock_true, L1_frequency);
        double l2_pseudo_range = gnss_satellites_.GetPseudoRangeECI(gnss_sat_id, sat_position_i, sat_clock_true, L2_frequency);
        auto l1_carrier_phase = gnss_satellites_.GetCarrierPhaseECI(gnss_sat_id, sat_position_i, sat_clock_true, L1_frequency);
        auto l2_carrier_phase = gnss_satellites_.GetCarrierPhaseECI(gnss_sat_id, sat_position_i, sat_clock_true, L2_frequency);

        double ionfree_range = (pow(L1_frequency/L2_frequency, 2.0)*l1_pseudo_range - l2_pseudo_range)/(pow(L1_frequency/L2_frequency, 2.0) - 1);
        double ionfree_phase = (pow(L1_frequency/L2_frequency, 2.0)*L1_lambda*(l1_carrier_phase.first + l1_carrier_phase.second) - L2_lambda*(l2_carrier_phase.first + l2_carrier_phase.second))/(pow(L1_frequency/L2_frequency, 2.0) - 1);

        gnss_true.can_see_satellites_sat_id.push_back(gnss_sat_id);
        gnss_true.gnss_satellites_position.push_back(gnss_position);
        gnss_true.gnss_clock.push_back(gnss_clock);
        gnss_true.L1_pseudo_range.push_back(l1_pseudo_range);
        gnss_true.L2_pseudo_range.push_back(l2_pseudo_range);
        gnss_true.L1_carrier_phase.push_back(l1_carrier_phase);
        gnss_true.L2_carrier_phase.push_back(l2_carrier_phase);

        gnss_true.ionfree_pseudo_range.push_back(ionfree_range);
        gnss_true.ionfree_carrier_phase.push_back(ionfree_phase);

        //観測誤差を混ぜる
        std::normal_distribution<> pseudo_range_noise(0.0, pseudo_sigma);
        std::normal_distribution<> carrier_phase_noise(0.0, carrier_sigma);

        //estimateの情報
        gnss_position = gnss_satellites_.GetSatellitePositionEci(i);
        gnss_clock = gnss_satellites_.GetSatelliteClock(i);

        l1_pseudo_range += pseudo_range_noise(mt); //pseudo_range_noise(mt);
        l2_pseudo_range += pseudo_range_noise(mt); //pseudo_range_noise(mt);

        l1_carrier_phase.first += carrier_phase_noise(mt)/L1_lambda; //carrier_phase_noise(mt)/L1_lambda;
        l2_carrier_phase.first += carrier_phase_noise(mt)/L2_lambda; //carrier_phase_noise(mt)/L2_lambda;

        ionfree_range = (pow(L1_frequency/L2_frequency, 2.0)*l1_pseudo_range - l2_pseudo_range)/(pow(L1_frequency/L2_frequency, 2.0) - 1);
        ionfree_phase = L2_lambda*(L1_frequency/L2_frequency* (l1_carrier_phase.first + l1_carrier_phase.second) - (l2_carrier_phase.first + l2_carrier_phase.second))/(pow(L1_frequency/L2_frequency, 2.0) - 1);

        gnss_observed_values.can_see_satellites_sat_id.push_back(gnss_sat_id);
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

void PBD_dgps::ProcessGnssObservation(GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, const int sat_id)
{
    int ind = 0;
    for(int i = 0;i < num_of_gnss_sastellites;++i){
        if(pre_observed_vector[sat_id].at(i) == true && now_observed_vector[sat_id].at(i) == false){
            first_L1_bias[sat_id].at(i) = 0.0;
            first_L2_bias[sat_id].at(i) = 0.0;
        }else if(pre_observed_vector[sat_id].at(i) == false && now_observed_vector[sat_id].at(i) == true){
            first_L1_bias[sat_id].at(i) = gnss_true.L1_carrier_phase.at(ind).second;
            first_L2_bias[sat_id].at(i) = gnss_true.L2_carrier_phase.at(ind).second;
        }
        if(now_observed_vector[sat_id].at(i)) ++ind;
    }

    for(int i = 0;i < gnss_true.can_see_satellites_sat_id.size();++i){
        int gnss_sat_id = gnss_true.can_see_satellites_sat_id.at(i);

        auto L1_observed = gnss_true.L1_carrier_phase.at(i);
        L1_observed.first += L1_observed.second - first_L1_bias[sat_id].at(gnss_sat_id);

        auto L2_observed = gnss_true.L2_carrier_phase.at(i);
        L2_observed.first += L2_observed.second - first_L2_bias[sat_id].at(gnss_sat_id);

        double ionfree_phase = (pow(L1_frequency/L2_frequency, 2.0)*L1_lambda*L1_observed.first - L2_lambda*L2_observed.first)/(pow(L1_frequency/L2_frequency, 2.0) - 1);
        gnss_true.ionfree_carrier_phase.push_back(ionfree_phase);
    }

    for(int i = 0;i < gnss_observed_values.can_see_satellites_sat_id.size();++i){
        int gnss_sat_id = gnss_observed_values.can_see_satellites_sat_id.at(i);

        auto L1_observed = gnss_observed_values.L1_carrier_phase.at(i);
        L1_observed.first += L1_observed.second - first_L1_bias[sat_id].at(gnss_sat_id);

        auto L2_observed = gnss_observed_values.L2_carrier_phase.at(i);
        L2_observed.first += L2_observed.second - first_L2_bias[sat_id].at(gnss_sat_id);

        double ionfree_phase = (pow(L1_frequency/L2_frequency, 2.0)*L1_lambda*L1_observed.first - L2_lambda*L2_observed.first)/(pow(L1_frequency/L2_frequency, 2.0) - 1);
        gnss_observed_values.ionfree_carrier_phase.push_back(ionfree_phase);
    }

    return;
}

void PBD_dgps::find_common_observed_gnss(const std::pair<int, int> sat_id_pair)
{
  // 初期化
  common_observed_gnss_sat_id.at(0).clear();
  common_observed_vector.at(0).assign(num_of_gnss_sastellites, false);
  const int main_sat_id = sat_id_pair.first;
  const int target_sat_id = sat_id_pair.second;
  int common_index = 0;
  for (int i = 0; i < now_observed_gnss_sat_id.at(main_sat_id).size(); ++i)
  {
    for (int j = 0; j < now_observed_gnss_sat_id.at(target_sat_id).size(); ++j)
    {
      int gnss_sat_id = now_observed_gnss_sat_id.at(main_sat_id).at(i);
      // どっかで複数衛星にも拡張
      if (gnss_sat_id == now_observed_gnss_sat_id.at(target_sat_id).at(j))
      {
        common_observed_gnss_sat_id[0].push_back(gnss_sat_id);
        common_observed_vector[0].at(gnss_sat_id) = true;
        //ここでいいのか不安
        common_index_dict.at(gnss_sat_id) = common_index;
        ++common_index;
        break;
      }
    }
  }
}

void PBD_dgps::calculate_phase_bias(GnssObservedValues gnss_observed_values, GnssObservedValues gnss_true, const int sat_id, Eigen::MatrixXd& pre_M, Eigen::VectorXd& bias)
{
  //観測する衛星同じだったら飛ばしていい
  if (check_vector_equal(pre_observed_gnss_sat_id.at(sat_id), now_observed_gnss_sat_id.at(sat_id))) {
    for (int i = 0; i < num_of_gnss_sastellites; ++i) now_observed_vector.at(sat_id).at(i) = false;
    return;
  }
  Eigen::VectorXd& bias_ = bias;
  int n = now_observed_gnss_sat_id.at(0).size();
  int n_pre = estimated_differential_bias.size();
  Eigen::VectorXd pre_bias = bias_;
  bias_.resize(n);

  vector<int> pre_index_from_sat_id(num_of_gnss_sastellites, -1);
  int index = 0;
  for (auto x : pre_observed_gnss_sat_id.at(sat_id)) {
    pre_index_from_sat_id.at(x) = index;
    ++index;
  }

  vector<int> now_index_from_sat_it(num_of_gnss_sastellites, -1);
  index = 0;
  for (auto x : now_observed_gnss_sat_id.at(sat_id)) {
    now_index_from_sat_it.at(x) = index;
    ++index;
  }

  int pre_index = 0; //index is not id
  int now_index = 0; //index is not id

  double geo_ure_sigma = sqrt(pow(L1_frequency / L2_frequency, 4.0) + 1.0) * pseudo_sigma / (pow(L1_frequency / L2_frequency, 2.0) - 1.0); //ionfree_pseudo_sigma

  int M_offset = (n + num_of_single_status) * sat_id;
  int pre_M_offset = (n_pre + num_of_single_status) * sat_id;
  // i = gnss_sat_id
  for (int i = 0; i < num_of_gnss_sastellites; ++i) {
    if (pre_observed_vector.at(sat_id).at(i) == false && now_observed_vector.at(sat_id).at(i) == false) continue;
    else if (pre_observed_vector.at(sat_id).at(i) == true && now_observed_vector.at(sat_id).at(i) == false) {
      if (pre_index != pre_index_from_sat_id.at(i)) {
        cout << "pre index something is wrong" << endl;
        abort();
      }
      ++pre_index;
    }
    else if (pre_observed_vector.at(sat_id).at(i) == false && now_observed_vector.at(sat_id).at(i) == true) {
      //前まではなかったので、対角成分だけどうにかすればいい
      bias_(now_index) = gnss_observed_values.ionfree_carrier_phase.at(now_index) - gnss_observed_values.ionfree_pseudo_range.at(now_index); // ここで入れてるバイアスは電離圏フリー結合のものになっている．っぽい

      //対角成分
      M(num_of_single_status + M_offset + now_index, num_of_single_status + M_offset + now_index) = pow(geo_ure_sigma, 2.0);

      if (now_index != now_index_from_sat_it.at(i)) {
        cout << "now index something is wrong 1" << endl;
        abort();
      }
      ++now_index;
    }
    else if (pre_observed_vector.at(sat_id).at(i) == true && now_observed_vector.at(sat_id).at(i) == true) {
      //バイアス成分
      bias_(now_index) = pre_bias(pre_index);

      //estimated_satusは引継ぎ
      M.block(0 + M_offset, now_index + num_of_single_status + M_offset, num_of_single_status, 1) = pre_M.block(0 + pre_M_offset, pre_index + num_of_single_status + pre_M_offset, num_of_single_status, 1);
      M.block(now_index + num_of_single_status + M_offset, 0 + M_offset, 1, num_of_single_status) = pre_M.block(pre_index + num_of_single_status + pre_M_offset, 0 + pre_M_offset, 1, num_of_single_status);

      //jは今存在するところしか走らない
      for (int j = 0; j < now_observed_gnss_sat_id.at(sat_id).size(); ++j) {
        if (j == now_index) break; //自分に追いついたのでbreak
        int j_now_gnss_sat_id = now_observed_gnss_sat_id.at(sat_id).at(j); //向こうの今のsat_id
        // いまいち何してんのかわからん<-MのNに関する成分をコピっている．関数化したい．
        if (pre_index_from_sat_id.at(j_now_gnss_sat_id) == -1) continue;
        int j_pre_index = pre_index_from_sat_id.at(j_now_gnss_sat_id);

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

      ++pre_index;   ++now_index;
    }
  }

  for (int i = 0; i < num_of_gnss_sastellites; ++i) {
    pre_observed_vector.at(sat_id).at(i) = now_observed_vector.at(sat_id).at(i);
    now_observed_vector.at(sat_id).at(i) = false;
  }
  pre_observed_gnss_sat_id.at(sat_id).clear();
  for (auto x : now_observed_gnss_sat_id.at(sat_id)) pre_observed_gnss_sat_id.at(sat_id).push_back(x);
}


void PBD_dgps::resize_Matrix(std::pair<GnssObservedValues, GnssObservedValues> gnss_observed_pair, std::pair<GnssObservedValues, GnssObservedValues> gnss_true_pair)
{
    // 共通衛星見つける
    find_common_observed_gnss(std::make_pair(0, 1));

    //バイアスの形状を変更
    Eigen::MatrixXd pre_M = M;
    int n_main = now_observed_gnss_sat_id.at(0).size();
    int n_main_pre = pre_observed_gnss_sat_id.at(0).size();
    int n_common = common_observed_gnss_sat_id.at(0).size();
    

    //M.resize(num_of_status + n_main + n_common, num_of_status + n_main + n_common);
    M = Eigen::MatrixXd::Zero(num_of_status + n_main + n_common, num_of_status + n_main + n_common); //0クリア
    //デフォルトのstatusはそのまま
    M.topLeftCorner(num_of_single_status, num_of_single_status) = pre_M.topLeftCorner(num_of_single_status, num_of_single_status);
    M.block(num_of_single_status + n_main, num_of_single_status + n_main, num_of_single_status, num_of_single_status) = pre_M.block(num_of_single_status + n_main_pre, num_of_single_status + n_main_pre, num_of_single_status, num_of_single_status);

    calculate_phase_bias(gnss_observed_pair.first, gnss_true_pair.first, 0, pre_M, estimated_bias);
    calculate_phase_bias(gnss_observed_pair.second, gnss_true_pair.second, 1, pre_M, estimated_target_bias);
    
    // differential_biasを求める．
    vector<int> target_index_from_sat_id(num_of_gnss_sastellites, -1);
    int index = 0;
    for (auto x : now_observed_gnss_sat_id.at(1)) {
      target_index_from_sat_id.at(x) = index;
      ++index;
    }

    int main_index = 0;
    int common_index = 0;
    int target_index = 0;
    estimated_differential_bias.resize(common_observed_gnss_sat_id.at(0).size());
    for (int i = 0; i < num_of_gnss_sastellites; ++i) {
      if (common_observed_vector.at(0).at(i))
      {
        main_index = main_index_dict.at(i);
        common_index = common_index_dict.at(i);
        target_index = target_index_from_sat_id.at(i);
        //バグってたらここが怪しい
        estimated_differential_bias(common_index) = estimated_target_bias(target_index) - estimated_bias(main_index);
      }
    }
    if (estimated_differential_bias.size() != common_observed_gnss_sat_id.at(0).size()) 
    {
      cout << "common index something wrong!" << endl;
      abort();
    }

    return;
}

bool PBD_dgps::CheckCanSeeSatellite(const libra::Vector<3> satellite_position, const libra::Vector<3> gnss_position) const
{
    double angle_rad = angle(satellite_position, gnss_position - satellite_position);
    if(angle_rad < M_PI/2.0 - mask_angle) return true;
    else return false;
}

// この関数がconstやからここでノイズを入れ込むことはできない？constはずそか．
double PBD_dgps::calculate_pseudo_range(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position, double gnss_clock) const // 本来はここのgnss_positionをIGSとかとIGUとかで分けて扱う必要がある
{
    double res = 0.0;
    res = calculate_geometric_range(sat_status, gnss_position);
    // clock offsetの分を追加
    res += sat_status(3) - gnss_clock; // 電離層はフリーにしている．
    // 観測ノイズ
    //std::normal_distribution<> pseudo_range_noise(0.0, pseudo_sigma);
    //res += pseudo_range_noise(mt);
    return res;
}

// 電離圏とかはまた考える．
double PBD_dgps::calculate_carrier_phase(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position, double gnss_clock, double integer_bias, double lambda) const
{
  double res = 0.0;
  res = calculate_geometric_range(sat_status, gnss_position);
  res += sat_status(3) - gnss_clock;
  // res -= lambda * integer_bias; // この形でできるようにするためにはbiasがちゃんと整数不定性の量になっている必要がある．
  res += integer_bias;
  // 観測ノイズ
  //std::normal_distribution<> carrier_phase_noise(0.0, carrier_sigma);
  // res += carrier_phase_noise(mt);

  return res;
}

double PBD_dgps::calculate_geometric_range(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position) const
{
    double res = 0.0;
    for(int i = 0;i < 3;++i){
        res += pow(sat_status(i) - gnss_position[i], 2.0);
    }
    res = sqrt(res);

    return res;
}

//使ってない，修正してこの形に持っていく
Eigen::VectorXd PBD_dgps::calculate_single_difference(const Eigen::VectorXd& main_observation, const Eigen::VectorXd& target_observation) const
{
  return main_observation - target_observation;
}


template <typename T> bool PBD_dgps::check_vector_equal(const vector<T>& a, const vector<T>& b)
{
    if(a.size() != b.size()) return false;
    for(int i = 0;i < a.size();++i){
        if(a.at(i) != b.at(i)) return false;
    }

    return true;
}