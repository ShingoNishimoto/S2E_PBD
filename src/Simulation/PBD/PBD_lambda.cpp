#include "PBD_const.h"
#include "PBD_Lambda.h"

PBD_Lambda::PBD_Lambda(Eigen::MatrixXd Qa, Eigen::VectorXd a_est)
{
  n = Qa.size();
  Z = Eigen::MatrixXd::Identity(n, n);
}

PBD_Lambda::~PBD_Lambda() {}


