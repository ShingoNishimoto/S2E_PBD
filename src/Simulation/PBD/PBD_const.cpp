#include "PBD_const.h"
#include <GnssSatellites.h>
#include <Environment/Global/PhysicalConstants.hpp>

// TODO: これも初期設定的なものはiniファイルに入れる．

const double L1_frequency = 1575.42; //[MHz]
const double L2_frequency = 1227.60;
const double L1_lambda = environment::speed_of_light_m_s*1e-6/L1_frequency; //[m]
const double L2_lambda = environment::speed_of_light_m_s*1e-6/L2_frequency; //[m]

const double mu_e_spice = 398600435436095.94;
const double mu_e = mu_e_spice;
// const double mu_e = environment::earth_gravitational_constant_m3_s2; //GM_E m^3/s^2
const double J2_const = 1.082616e-3; // 無次元 重力J2項
const double earth_radius_spice_mean = 6371000.3852496156;
// const double Earth_Radius = earth_radius_spice_mean;
const double Earth_Radius = environment::earth_equatorial_radius_m; //m
