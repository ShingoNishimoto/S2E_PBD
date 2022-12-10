#include "PBD_const.h"
#include <GnssSatellites.h>
#include <Environment/Global/PhysicalConstants.hpp>

// TODO: これも初期設定的なものはiniファイルに入れる．

const double L1_frequency = 1575.42; //[MHz]
const double L2_frequency = 1227.60;
const double L1_lambda = environment::speed_of_light_m_s*1e-6/L1_frequency; //[m]
const double L2_lambda = environment::speed_of_light_m_s*1e-6/L2_frequency; //[m]

// A-priori standard deviation
const double sigma_r_ini     = 100;   //[m]
const double sigma_v_ini     = 1.0;  //[m/s]
const double sigma_acc_r_ini =  75*1e2; //[nm/s^2]
const double sigma_acc_t_ini = 100*1e2; //[nm/s^2]
const double sigma_acc_n_ini =  50*1e2;  //[nm/s^2]
const double sigma_cdt_ini   = 100;    //[m]
const double sigma_N_ini     = 10;   //[cycle] これももう少し現実にそった値にする．
// observation noise
// const double pseudo_sigma    = 1.0 / 3.0; //[m]
// const double carrier_sigma   = 5.0 * 1e-3 / 3.0; //[m]
// for SILVIA
const double pseudo_sigma    = 2.0 / 3; //[m]
const double carrier_sigma   = 4.0 * 1e-3 / 3; //[m]
// process noise
const double sigma_r_process = 0.001;    //[m]
const double sigma_v_process = 5 * 1e-4;   //[m/s]
const double sigma_acc_r_process = 50*1e2; // 7500;   //[nm/s^2]
const double sigma_acc_t_process = 80*1e2; // 1000;   //[nm/s^2]
const double sigma_acc_n_process = 40*1e2; // 500;    //[nm/s^2]
// DiGiTaLでは500mになっている．実際の受信機使っているからか？
const double sigma_cdt_process   = 1.0; // 0.25;    //[m] <- これもホンマはホワイトノイズとランダムウォークに分ける必要がある．ドリフトと，バイアス．

const double sigma_N_process     = 0.1;    //[cycle]
// clock noise model parameters of receiver
const double clock_sigma     = 1.0; // 0.25; //[m] 0.1 で<- 0.1nsくらいになる．これは今white noiseになっている．

const double mu_e_spice = 398600435436095.94;
const double mu_e = mu_e_spice;
// const double mu_e = environment::earth_gravitational_constant_m3_s2; //GM_E m^3/s^2
const double J2_const = 1.082616e-3; // 無次元 重力J2項
const double earth_radius_spice_mean = 6371000.3852496156;
// const double Earth_Radius = earth_radius_spice_mean;
const double Earth_Radius = environment::earth_equatorial_radius_m; //m

const double tau_a = 900;
const double tau_cdt = 100;

// for Adaptive Kalman Filter
const double alpha = 0.3; // forgetting factor
