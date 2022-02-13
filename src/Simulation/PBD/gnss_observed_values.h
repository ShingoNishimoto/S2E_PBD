#ifndef GNSS_OBSERVED_VALUES_H_
#define GNSS_OBSERVED_VALUES_H_

#include <vector>
#include <utility>
#include "Vector.hpp"

struct GnssObservedValues
{
    std::vector<int> observable_gnss_sat_id;
    std::vector<libra::Vector<3>> gnss_satellites_position;
    std::vector<double> gnss_clock;

    std::vector<double> L1_pseudo_range; //[m]
    std::vector<double> L2_pseudo_range; //[m]
    std::vector<std::pair<double, double>> L1_carrier_phase; //[位相, bias(整数, dtype = double)]
    std::vector<std::pair<double, double>> L2_carrier_phase; //[位相, bias(整数, dtype = double)]

    std::vector<double> ionfree_pseudo_range;
    std::vector<double> ionfree_carrier_phase; //今入っているのは加工済み

    bool check_normal();
};

#endif
