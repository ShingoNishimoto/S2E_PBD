#ifndef __PBD_pod_H__
#define __PBD_pod_H__

#include "Dense" // Eigen

#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include "Vector.hpp"
#include "SimTime.h"
#include "GnssSatellites.h"
#include "./Orbit/Orbit.h"

#include "gnss_observed_values.h"

class PBD_pod
{       
    public:
        PBD_pod(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Orbit& orbit_);
        ~PBD_pod();
        void Update(const SimTime& sim_time_, const GnssSatellites& gnss_satellites_, const Orbit& orbit_);
        void OrbitPropagation();
        void GetGnssPositionObservation(const GnssSatellites& gnss_satellites_, const Orbit& orbit_, GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true);
        void ProcessGnssObservation(GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true);
        void resize_Matrix(GnssObservedValues& gnss_observed_values, GnssObservedValues& gnss_true);
        void KalmanFilter(const GnssObservedValues& gnss_observed_values, const GnssObservedValues& gnss_true);

        //仮
        //アンテナの中心の向きが、常に反地球方向を向いているとして、適当にマスク角を取って、その中にいるとする
        bool CheckCanSeeSatellite(const libra::Vector<3> satellite_position, const libra::Vector<3> gnss_position) const;

        std::ofstream ofs;

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    private:
        //[x:m, y:m, z:m, s:[m]]
        Eigen::Vector4d estimated_status;
        //[vx[m/s], vy[m/s], vz[m/s]]
        Eigen::Vector3d estimated_velocity;
        //[ax[m/s^2], ay[m/s^2], az[m/s^2]]
        Eigen::Vector3d estimated_acc;

        Eigen::VectorXd estimated_bias;

        vector<bool> pre_observed_vector;
        vector<bool> now_observed_vector;
        vector<int> pre_observed_sat_id;
        vector<int> now_observed_sat_id;

        vector<double> first_L1_bias;
        vector<double> first_L2_bias;
        int num_of_sastellites;

        const double Cd = 2.928e-14;

        Eigen::MatrixXd M;

        const int num_of_status = 7;

        double step_time;
        double observe_step_time = 10.0;
        double log_step_time = 1.0;
        Eigen::Vector3d position_differential(const Eigen::Vector3d& velocity) const;
        Eigen::Vector3d velocity_differential(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const;
        Eigen::MatrixXd update_M_matrix(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity);
        Eigen::MatrixXd calculate_A_matrix(const Eigen::Vector3d& position, const Eigen::Vector3d& velocity) const;
        Eigen::MatrixXd calculate_Q_matrix(const int n);

        //clockの真値[m]
        double sat_clock_true;

        //マスク角 [rad]
        const double mask_angle = 10.0/180.0*M_PI;

        std::random_device seed_gen;
        std::mt19937 mt;

        double pseudo_range_calculator(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position, double gnss_clock) const;
        double geo_range_calculator(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position) const;
        double carrier_phase_calculator(const Eigen::Vector4d& sat_status, libra::Vector<3> gnss_position, double gnss_clock, double integer_bias, double lambda_narrow) const;

        template <typename T> bool check_vector_equal(const vector<T>& a, const vector<T>& b);
};
#endif