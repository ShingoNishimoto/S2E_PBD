#include "PBD_GeoPotential.hpp"
#include <InitDisturbance.hpp>

PBD_GeoPotential::PBD_GeoPotential(const int degree, const LocalEnvironment& local_env, std::string file_path) : degree_(degree), local_env_(local_env)
{
  GeoPotential geo_potential(degree, file_path);
  geo_potential_ = &geo_potential;
  geo_potential_->IsCalcEnabled = false;
  geo_potential_->IsLogEnabled = false; // ログも多分要らない．
}

PBD_GeoPotential::~PBD_GeoPotential() {}

Vector<3> PBD_GeoPotential::CalcAccelerationECI(const Vector<3> &position_eci)
{
  Matrix<3, 3> trans_eci2ecef_ = local_env_.GetCelesInfo().GetGlobalInfo().GetEarthRotation().GetDCMJ2000toXCXF();
  const Vector<3> position_ecef = trans_eci2ecef_ * position_eci;

  geo_potential_->CalcAccelerationECEF(position_ecef);
  // acc_ecef_がprivateの変数やから微妙．このままではアクセスできない．継承しても無理なので同じコードを爆誕させるしかないか．一旦，coreを修正する．
  const Vector<3> acc_ecef_ =  geo_potential_->GetAccECEF();
  Matrix<3, 3> trans_ecef2eci_ = transpose(trans_eci2ecef_);

  Vector<3> acc_eci_ = trans_ecef2eci_ * acc_ecef_;
  return acc_eci_;
}

