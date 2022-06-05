#include "PBD_const.h"
#include "PBD_Lambda.h"

PBD_Lambda::PBD_Lambda(Eigen::MatrixXd Q_a, Eigen::VectorXd a_est): a(a_est)
{
  n = Q_a.size();
  Z = Eigen::MatrixXd::Identity(n, n);
  // Ç†Ç¡ÇƒÇ¢ÇÈÇÃÇ©ÅH
  L = Eigen::MatrixXd::Zero(n, n);
  L.triangularView<Eigen::UnitLower>() = Eigen::MatrixXd::Identity(n, n);
}

PBD_Lambda::~PBD_Lambda() {}


void PBD_Lambda::IntegerGaussTransformation()
{
  int mu;
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; i < n; ++i)
    {
      mu = L(i, j);
      if (mu != 0)
      {
        L.block(i, j, n - i, 1) = L.block(i, j, n - i, 1) - mu * L.block(i, i, n - i, 1);
        Z.block(0, j, n, 1) = Z.block(0, j, n, 1) - mu * Z.block(0, i, n, 1);
        a(j) = a(j) - mu * a(i);
      }
    }
  }

}
