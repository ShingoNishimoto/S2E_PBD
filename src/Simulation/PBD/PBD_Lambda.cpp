#include "PBD_const.h"
#include "PBD_Lambda.h"
#include <iostream>

#define LAMBDA_DEBUG

PBD_Lambda::PBD_Lambda(Eigen::MatrixXd Q_a, Eigen::VectorXd a_est): a(a_est)
{
  n = Q_a.rows();
  Z = Eigen::MatrixXd::Identity(n, n);
  Eigen::LDLT<Eigen::MatrixXd> LDLTOfQ_a(Q_a.transpose());
  L = LDLTOfQ_a.matrixL();
  D = LDLTOfQ_a.vectorD();
  // L.triangularView<Eigen::UnitLower>() = Eigen::MatrixXd::Identity(n, n);
}

PBD_Lambda::~PBD_Lambda() {}

void PBD_Lambda::Solve(void)
{
  IntegerGaussTransformation();
  Permute();
  return;
}

void PBD_Lambda::IntegerGaussTransformation(void)
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
#ifdef LAMBDA_DEBUG
  std::cout << "L" << L << std::endl;
  std::cout << "Z" << Z << std::endl;
  std::cout << "a" << a << std::endl;
  std::cout << "LZ" << L*Z << std::endl;
#endif // LAMBDA_DEBUG
}

void PBD_Lambda::Permute(void)
{
  for (int k = 0; k < n - 1; ++k)
  {
    double delta = D(k) + L(k + 1, k) * L(k + 1, k) * D(k + 1);
    double ita = D(k) / delta;
    double lambda = D(k + 1) * L(k + 1, k) / delta;
    D(k) = ita * D(k + 1);
    D(k) = delta;
    Eigen::Matrix2d temp;
    temp << -L(k + 1, k), 1,
             ita, lambda;
    if (k > 1) L.block(k, 0, 2, k - 1) = temp * L.block(k, 0, 2, k - 1);
    L(k + 1, k) = lambda;
    // L swap operation
    Z.col(k).swap(Z.col(k + 1));
    double a_k = a(k);
    a(k+1) = a(k + 1);
    a(k) = a_k;
  }
#ifdef LAMBDA_DEBUG
  std::cout << "D" << D << std::endl;
  std::cout << "Z" << Z << std::endl;
  std::cout << "a" << a << std::endl;
#endif // LAMBDA_DEBUG
}

