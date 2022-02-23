#include "PBD_const.h"
#include "GnssSatellites.h"

const double L1_frequency = 1575.42; //[MHz]
const double L2_frequency = 1227.60;
const double L1_lambda = speed_of_light*1e-6/L1_frequency; //[m]
const double L2_lambda = speed_of_light*1e-6/L2_frequency; //[m]

// A-priori standard deviation
const double sigma_r_ini     = sqrt(10);   //[m]
const double sigma_v_ini     = 1.0;  //[m/s]
const double sigma_acc_r_ini = 1000; //[nm/s^2]
const double sigma_acc_t_ini = 2000; //[nm/s^2]
const double sigma_acc_n_ini = 750;  //[nm/s^2]
const double sigma_cdt_ini   = 1;    //[m]
const double sigma_N_ini     = 1;    //[m]
// observation noise
const double pseudo_sigma    = 1; //[m]
const double carrier_sigma   = 1e-2; //[m]
// process noise
const double sigma_r_process = 0.1;    //[m]
const double sigma_v_process = 1e-1;    //[m/s]
// 上の二つは入れていいのかを要検討
const double sigma_acc_r_process = 7500;    //[nm/s^2]
const double sigma_acc_t_process = 1000;    //[nm/s^2]
const double sigma_acc_n_process = 500;    //[nm/s^2]
const double sigma_cdt_process   = 0.1;    //[m] <- これもホンマはホワイトノイズとランダムウォークに分ける必要がある．ドリフトと，バイアス．　
const double sigma_N_process     = 0.1;    //[m]
// clock noise model parameters
const double clock_sigma     = 0.1;  //[m] 0.1 で<- 0.1nsくらいになる．これはrandom walk

const double mu_const = 3.986004418e14; //GM_E m^3/s^2
const double J2_const = 1.082636e-3; //無次元 重力J2項
const double Earth_Radius = 6378136.6; //m

const double tau_a = 400;
const double tau_cdt = 100;
