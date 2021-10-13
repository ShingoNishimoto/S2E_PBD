#include "PBD_const.h"
#include "GnssSatellites.h"

const double L1_frequency = 1575.42; //[MHz]
const double L2_frequency = 1227.60;
const double L1_lambda = speed_of_light*1e-6/L1_frequency; //[m];
const double L2_lambda = speed_of_light*1e-6/L2_frequency; //[m]
const double pseudo_sigma = 1.0;
const double carrier_sigma = 1e-3;
const double clock_sigma = 1.0;

const double mu_const = 3.986004418e14; //GM_E m^3/s^2
const double J2_const = 1.082636e-3; //無次元 重力J2項
const double Earth_Radius = 6378136.6; //m
