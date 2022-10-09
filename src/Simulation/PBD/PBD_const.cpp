#include "PBD_const.h"
#include <GnssSatellites.h>
#include <Environment/Global/PhysicalConstants.hpp>

// TODO: これも初期設定的なものはiniファイルに入れる．

const double L1_frequency = 1575.42; //[MHz]
const double L2_frequency = 1227.60;
const double L1_lambda = environment::speed_of_light_m_s*1e-6/L1_frequency; //[m]
const double L2_lambda = environment::speed_of_light_m_s*1e-6/L2_frequency; //[m]

// A-priori standard deviation
const double sigma_r_ini     = 1000;   //[m]
const double sigma_v_ini     = 1.0;  //[m/s]
const double sigma_acc_r_ini = 1000; // 150; //[nm/s^2]
const double sigma_acc_t_ini = 2000; // 300; //[nm/s^2]
const double sigma_acc_n_ini = 750; // 200;  //[nm/s^2]
const double sigma_cdt_ini   = 100;    //[m]
const double sigma_N_ini     = 10.0;   //[cycle]
// observation noise
const double pseudo_sigma    = 0.25; //[m]
const double carrier_sigma   = 5*1e-3; //[m]
// process noise
const double sigma_r_process = 0.001;    //[m]
const double sigma_v_process = 5 * 1e-4;   //[m/s]
const double sigma_acc_r_process = 500; // 7500;   //[nm/s^2]
const double sigma_acc_t_process = 200; // 1000;   //[nm/s^2]
const double sigma_acc_n_process = 100; // 500;    //[nm/s^2]
// DiGiTaLでは500mになっている．実際の受信機使っているからか？
const double sigma_cdt_process   = 5;    //[m] <- これもホンマはホワイトノイズとランダムウォークに分ける必要がある．ドリフトと，バイアス．

const double sigma_N_process     = 0.10;    //[cycle]
// clock noise model parameters of receiver
const double clock_sigma     = 5;  //[m] 0.1 で<- 0.1nsくらいになる．これは今white noiseになっている．

const double mu_e = 3.986004418e14; //GM_E m^3/s^2
const double J2_const = 1.082630e-3; // 無次元 重力J2項
const double Earth_Radius = environment::earth_equatorial_radius_m; //m

const double tau_a = 900;
const double tau_cdt = 100;

// for Adaptive Kalman Filter
const double alpha = 0.3; // forgetting factor
