#include "gnss_observed_values.h"

bool GnssObservedValues::check_normal()
{
    unsigned int n = can_see_satellites_sat_id.size();
    if(n != gnss_satellites_position.size()) return false;
    if(n != gnss_clock.size()) return false;
    if(n != L1_pseudo_range.size()) return false;
    if(n != L2_pseudo_range.size()) return false;
    if(n != L1_carrier_phase.size()) return false;
    if(n != L2_carrier_phase.size()) return false;
    if(n != ionfree_pseudo_range.size()) return false;
    //if(n != ionfree_carrier_phase.size()) return false;

    return true;
}