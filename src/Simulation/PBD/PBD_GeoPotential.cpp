#include "PBD_GeoPotential.hpp"
#include <InitDisturbance.hpp>

PBD_GeoPotential::PBD_GeoPotential(const int degree, std::string file_path)
{
  geop_ = new GeoPotential(degree, file_path);
  geop_->IsCalcEnabled = false;
  geop_->IsLogEnabled = false; // ログも多分要らない．
}

PBD_GeoPotential::~PBD_GeoPotential() {}

Vector<3> PBD_GeoPotential::CalcAccelerationECI(const Vector<3> &position_eci, Matrix<3, 3> trans_eci_to_ecef)
{
  const Vector<3> position_ecef = trans_eci_to_ecef * position_eci;

  geop_->CalcAccelerationECEF(position_ecef);
  // acc_ecef_がprivateの変数やから微妙．このままではアクセスできない．継承しても無理なので同じコードを爆誕させるしかないか．一旦，coreを修正する．
  const Vector<3> acc_ecef_ =  geop_->GetAccECEF();
  Matrix<3, 3> trans_ecef2eci_ = transpose(trans_eci_to_ecef);

  Vector<3> acc_eci_ = trans_ecef2eci_ * acc_ecef_;
  return acc_eci_;
}

