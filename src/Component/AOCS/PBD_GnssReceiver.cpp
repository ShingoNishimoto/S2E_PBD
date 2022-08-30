#include "PBD_GnssReceiver.hpp"
#include "../../Simulation/Spacecraft/PBD_Components.h"

void PBD_GnssReceiver::MainRoutine(int count)
{
  UNUSED(count);

  Vector<3> pos_true_eci_ = dynamics_->GetOrbit().GetSatPosition_i();
  Quaternion q_i2b = dynamics_->GetQuaternion_i2b();

  // 特定のGNSS衛星に対してみているのではなくて，反地球方向を見ているかの簡単な確認．
  CheckAntenna(pos_true_eci_, q_i2b);

  if (is_gnss_sats_visible_ == 1) {  // Antenna of GNSS-R can detect GNSS signal
    position_ecef_ = dynamics_->GetOrbit().GetSatPosition_ecef();
    position_llh_ = dynamics_->GetOrbit().GetLatLonAlt();
    velocity_ecef_ = dynamics_->GetOrbit().GetSatVelocity_ecef();
    // AddNoise(pos_true_eci_, position_ecef_); noiseをこっちでは追加しない．

    utc_ = simtime_->GetCurrentUTC();
    ConvertJulianDayToGPSTime(simtime_->GetCurrentJd());
  } else {
    // position information will not be updated in this case
    // only time information will be updated in this case (according to the
    // receiver's internal clock)
    utc_ = simtime_->GetCurrentUTC();
    ConvertJulianDayToGPSTime(simtime_->GetCurrentJd());
  }
}
