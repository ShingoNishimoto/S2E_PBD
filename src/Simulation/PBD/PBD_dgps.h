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

#include "gnss_observed_values.h"

class PBD_dgps
{       
    public:
        PBD_dgps(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Orbit& main_orbit, const Orbit& target_orbit);
        ~PBD_dgps();
        void Update(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_);//, const Orbit& main_orbit, const Orbit& target_orbit);
        void OrbitPropagation();
        void GetGnssPositionObservation(const GnssSatellites& gnss_satellites_, const Orbit& orbit_, const int sat_id, GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, double sat_clock_true);
        void ProcessGnssObservation(GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, const int sat_id);
        void resize_Matrix(std::pair<GnssObservedValues, GnssObservedValues> gnss_observed_pair, std::pair<GnssObservedValues, GnssObservedValues> gnss_true_pair);
        void calculate_phase_bias(GnssObservedValues gnss_observed_values, GnssObservedValues gnss_true, const int sat_id, Eigen::MatrixXd& pre_M, Eigen::VectorXd& bias);
        //void calculate_difference_observation(GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true, const int sat_id, Eigen::MatrixXd& pre_M);
        void KalmanFilter(const GnssObservedValues& gnss_observed_main, const GnssObservedValues& gnss_observed_target);

        //仮
        //アンテナの中心の向きが、常に反地球方向を向いているとして、適当にマスク角を取って、その中にいるとする
        bool CheckCanSeeSatellite(const libra::Vector<3> satellite_position, const libra::Vector<3> gnss_position) const;

        std::ofstream ofs;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:
        // main satellite
        const Orbit& main_orbit_;
        // target satellite
        const Orbit& target_orbit_;
        //[x:m, y:m, z:m, s:[m]]
        Eigen::Vector4d estimated_status;
        //[vx[m/s], vy[m/s], vz[m/s]]
        Eigen::Vector3d estimated_velocity;
        //[ax[mm/s^2], ay[mm/s^2], az[mm/s^2]]
        Eigen::Vector3d estimated_acc; // ここの加速度は外乱とかのその他加速度

        Eigen::VectorXd estimated_bias;

        std::map<const int, int> main_index_dict;

        // ここもsingleとdoubleで分けて考えれるような枠組みにしたい
        //[dx[m], dy[m], dz[m], ds[m]]
        Eigen::Vector4d estimated_differential_status;
        //[dvx[m/s], dvy[m/s], dvz[m/s]] <- mm/sの方がよさそう．
        Eigen::Vector3d estimated_differential_velocity;
        //[dax[mm/s^2], day[mm/s^2], daz[mm/s^2]]
        Eigen::Vector3d estimated_differential_acc;
        // 辞書が欲しい commonとmainをつなげるために
        Eigen::VectorXd estimated_target_bias;
        Eigen::VectorXd estimated_differential_bias;

        std::map<const int, int> common_index_dict;
        
        //初期化をもっとスマートにできるように考える
        vector<vector<bool>> pre_observed_vector{ {},{} }; // 名前分かりにくいので変えたい
        vector<vector<bool>> now_observed_vector{ {},{} };
        vector<vector<bool>> common_observed_vector{ {},{} };
        vector<vector<int>> pre_observed_gnss_sat_id{ {},{} };
        vector<vector<int>> now_observed_gnss_sat_id{ {},{} };
        vector<vector<int>> common_observed_gnss_sat_id{ {},{} };

        vector<vector<double>> first_L1_bias{ {},{} };
        vector<vector<double>> first_L2_bias{ {},{} };
        int num_of_gnss_sastellites;

        const double Cd = 2.928e-14; // 高度に応じて変更したいが，高度変化ないから一定でいいか．

        Eigen::MatrixXd M;

        const int num_of_single_status = 10; //+3+ns
        const int num_of_status = 20; //+3　+ns

        double step_time;
        double observe_step_time = 10.0;
        double log_step_time = 1.0;
        vector<Eigen::Vector3d> RK4(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration);
        Eigen::Vector3d position_differential(const Eigen::Vector3d& velocity) const;
        Eigen::Vector3d velocity_differential(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, Eigen::Vector3d& acceleration) const;
        // for differential
        Eigen::MatrixXd update_M_matrix(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration, const Eigen::Vector3d& position_difference, const Eigen::Vector3d& velocity_difference, const Eigen::Vector3d& acceleration_difference);
        Eigen::MatrixXd calculate_A_matrix(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration) const;
        // for differential
        Eigen::MatrixXd calculate_A_matrix(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity, const Eigen::Vector3d& acceleration, const Eigen::Vector3d& position_difference, const Eigen::Vector3d& velocity_difference, const Eigen::Vector3d& acceleration_difference) const;
        Eigen::MatrixXd calculate_Q_matrix(const int n_main, const int n_common);
        Eigen::MatrixXd calculate_Phi_a(const double dt);
        void find_common_observed_gnss(const std::pair<int, int> sat_id_pair);


        //clockの真値[m]
        double sat_main_clock_true;
        double sat_target_clock_true;

        //マスク角 [rad]
        const double mask_angle = 10.0/180.0*M_PI;

        std::random_device seed_gen;
        std::mt19937 mt;

        double calculate_pseudo_range(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position, double gnss_clock) const;
        double calculate_carrier_phase(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position, double gnss_clock, double integer_bias, double lambda) const;
        double calculate_geometric_range(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position) const;
        Eigen::VectorXd calculate_single_difference(const Eigen::VectorXd& main_observation, const Eigen::VectorXd& target_observation) const;
        // 一旦singleだけにする
        //double calculate_double_difference(const Eigen::VectorXd& main_observation, const Eigen::VectorXd& target_observation) const;

        template <typename T> bool check_vector_equal(const vector<T>& a, const vector<T>& b);
};
#endif