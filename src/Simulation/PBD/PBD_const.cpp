#include "PBD_const.h"
#include <GnssSatellites.h>
#include <Environment/Global/PhysicalConstants.hpp>

const double L1_frequency = 1575.42; //[MHz]
const double L2_frequency = 1227.60;
const double L1_lambda = environment::speed_of_light_m_s*1e-6/L1_frequency; //[m]
const double L2_lambda = environment::speed_of_light_m_s*1e-6/L2_frequency; //[m]

// A-priori standard deviation
const double sigma_r_ini     = 10;   //[m]
const double sigma_v_ini     = 1.0;  //[m/s]
const double sigma_acc_r_ini = 150; // 10000; //[nm/s^2]
const double sigma_acc_t_ini = 300; // 20000; //[nm/s^2]
const double sigma_acc_n_ini = 200; // 15000;  //[nm/s^2]
const double sigma_cdt_ini   = 100;    //[m]
const double sigma_N_ini     = 1;    //[m]
// observation noise
const double pseudo_sigma    = 1; //[m]
const double carrier_sigma   = 1e-3; //[m]
// process noise
const double sigma_r_process = 0;    //[m]
const double sigma_v_process = 0;    //[m/s] 
const double sigma_acc_r_process = 15; // 4000;    //[nm/s^2]
const double sigma_acc_t_process = 30; // 8000;   //[nm/s^2]
const double sigma_acc_n_process = 20; // 5000;    //[nm/s^2]
// DiGiTaLでは500mになっている．実際の受信機使っているからか？
const double sigma_cdt_process   = 5;    //[m] <- これもホンマはホワイトノイズとランダムウォークに分ける必要がある．ドリフトと，バイアス．

const double sigma_N_process     = 0.1;    //[m]
// clock noise model parameters of receiver
const double clock_sigma     = 10;  //[m] 0.1 で<- 0.1nsくらいになる．これはrandom walk

const double mu_const = 3.986004418e14; //GM_E m^3/s^2
const double J2_const = 1.082636e-3; //無次元 重力J2項
const double Earth_Radius = 6378136.6; //m

const double tau_a = 400;
const double tau_cdt = 100;
