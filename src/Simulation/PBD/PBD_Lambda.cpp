#include "PBD_const.h"
#include "PBD_Lambda.h"
#include <iostream>
#include <cassert>
#define LAMBDA_DEBUG

PBD_Lambda::PBD_Lambda(Eigen::MatrixXd P_a, Eigen::VectorXd a_est, const int n_cands): a_hat(a_est), n_cands(n_cands)
{
  n = P_a.rows();
  Z = Eigen::MatrixXd::Identity(n, n);
  Eigen::LDLT<Eigen::MatrixXd> LDLTOfP_a(P_a.transpose());
  L = LDLTOfP_a.matrixL();
  D = LDLTOfP_a.vectorD();
  z = Eigen::MatrixXd::Zero(n, n_cands);
  sq_norm = Eigen::VectorXd::Zero(n_cands);
#ifdef LAMBDA_DEBUG
  std::cout << "P_a" << P_a << std::endl;
  std::cout << "a_est" << a_est << std::endl;
#endif // LAMBDA_DEBUG
  // a ?
  // L.triangularView<Eigen::UnitLower>() = Eigen::MatrixXd::Identity(n, n);
}

PBD_Lambda::~PBD_Lambda() {}

void PBD_Lambda::Solve(void)
{
  Reduction();
  int n_candidates = 2;
  Eigen::VectorXd z_hat = Z.transpose() * a_hat;
  IntegerLeastSquares(n_candidates, z_hat);
#ifdef LAMBDA_DEBUG
  std::cout << "D" << D << std::endl;
  std::cout << "Z" << Z << std::endl;
  std::cout << "L" << L << std::endl;
  std::cout << "a_hat" << a_hat << std::endl;
  std::cout << "a" << a << std::endl;
#endif // LAMBDA_DEBUG
  return;
}

void PBD_Lambda::IntegerGaussTransformation(const int i, const int j)
{
  int mu;
  mu = (int)round(L(i, j));
  if (mu != 0)
  {
    // ‚±‚ê‚Í¡SD‚É‘Î‚µ‚Ä‚â‚Á‚Ä‚¢‚é‚©‚ç–³‘ŠŠÖ‚Å‚±‚±‚É“ü‚ç‚È‚¢DQ‚ðŒ©‚éŠ´‚¶–³‘ŠŠÖ‚Å‚Í‚È‚³‚»‚¤‚â‚¯‚ÇC‚Ù‚Ú–³‘ŠŠÖ‚È‚Ì‚©D
    L.block(i, j, n - i, 1) = L.block(i, j, n - i, 1) - mu * L.block(i, i, n - i, 1);
    Z.block(0, j, n, 1) = Z.block(0, j, n, 1) - mu * Z.block(0, i, n, 1);
    a_hat(j) = a_hat(j) - mu * a_hat(i);
  }
#ifdef LAMBDA_DEBUG
  /*
  std::cout << "L" << L << std::endl;
  std::cout << "Z" << Z << std::endl;
  std::cout << "a" << a_hat << std::endl;
  std::cout << "LZ" << L*Z << std::endl;
  */
#endif // LAMBDA_DEBUG
}

// k‚Ì‹–‚³‚ê‚é”ÍˆÍ 1 < k < n-3
void PBD_Lambda::Permute(const int k, const double delta)
{
  double ita = D(k) / delta;
  double lambda = D(k + 1) * L(k + 1, k) / delta;
  D(k) = ita * D(k + 1);
  D(k) = delta;
  Eigen::Matrix2d temp;
  temp << -L(k + 1, k), 1,
            ita, lambda;
  if (k > 0) L.block(k, 0, 2, k) = temp * L.block(k, 0, 2, k);
  L(k + 1, k) = lambda;
  Eigen::MatrixXd L_temp = L.block(k + 2, k, n - (k + 2), 1);
  L.block(k + 2, k, n - (k + 2), 1) = L.block(k + 2, k + 1, n - (k + 2), 1);
  L.block(k + 2, k + 1, n - (k + 2), 1) = L_temp;
  Z.col(k).swap(Z.col(k + 1));
  double a_k = a_hat(k);
  a_hat(k+1) = a_hat(k + 1);
  a_hat(k) = a_k;

#ifdef LAMBDA_DEBUG
  // std::cout << "D" << D << std::endl;
  // std::cout << "Z" << Z << std::endl;
  // std::cout << "a" << a_hat << std::endl;
#endif // LAMBDA_DEBUG
}

void PBD_Lambda::Reduction(void)
{
  int k, l;
  k = n - 2;
  l = n - 2;
  Eigen::VectorXd D_bar = Eigen::VectorXd::Zero(n);
  while (k >= 0)
  {
    assert(k < n - 1);
    if (k <= l)
    {
      for (int i = k + 1; i < n; i++)
      {
        IntegerGaussTransformation(i, k);
      }
    }
    D_bar(k + 1) = D(k) + pow(L(k + 1, k), 2) * D(k + 1);
    if (D_bar(k + 1) < D(k + 1))
    {
      Permute(k, D_bar(k + 1));
      l = k;
      k = n - 2;
    }
    else
    {
      k--;
    }
  }
}

void PBD_Lambda::IntegerLeastSquares(const int n_cands, Eigen::VectorXd z_hat)
{
  assert(z_hat.size() > 0);

  double chi2 = 1.0 * 10e+18;
  Eigen::VectorXd dist = Eigen::VectorXd::Zero(n);
  int search_finished = 0;
  int count = -1; // the number of candidates
  
  Eigen::VectorXd z_cond = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd z_round = Eigen::VectorXd::Zero(n);
  Eigen::VectorXi step = Eigen::VectorXi::Zero(n);

  z_cond(n - 1) = z_hat(n - 1);
  z_round(n - 1) = round(z_cond(n - 1));

  double diff = z_cond(n - 1) - z_round(n - 1);
  step(n - 1) = (diff >= 0) ? 1 : -1;

  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
  int k = n - 1;
  int i_max = n_cands - 1; // Initially, maximum F(z) is at n_cands

  while (!search_finished)
  {
    double dist_new = dist(k) + pow(diff, 2) / D(k);
    if (dist_new < chi2)
    {
      if (k != 0)
      {
        k--;
        dist(k) = dist_new;
        S.block(k, 0, 1, k + 1) = S.block(k + 1, 0, 1, k + 1) + (z_round(k + 1) - z_cond(k + 1)) * L.block(k + 1, 0, 1, k + 1);

        z_cond(k) = z_hat(k) + S(k, k);
        z_round(k) = round(z_cond(k));
        diff = z_cond(k) - z_round(k);
        step(k) = (diff >= 0) ? 1 : -1;
      }
      else
      {
        // store the found candidate and try nect valid integer
        if (count < n_cands - 2)
        {
          count++;
          z.col(count) = z_round;
          sq_norm(count) = dist_new;
        }
        else
        {
          z.col(i_max) = z_round;
          sq_norm(i_max) = dist_new;

          i_max = Argmax(sq_norm);
          chi2 = sq_norm(i_max);
        }
        z_round(0) += step(0); // next valid integer
        diff = z_cond(0) - z_round(0);
        int step_sign = (step(0) > 0) ? 1 : (step(0) == 0) ? 0 : -1;
        step(0) = -step(0) - step_sign;
      }
    }
    else
    {
      // exit or move up
      if (k == n - 1) search_finished = 1;
      else
      {
        k++;
        z_round(k) += step(k); // next valid integer
        diff = z_cond(k) - z_round(k);
        int step_sign = (step(k) > 0) ? 1 : (step(k) == 0) ? 0 : -1;
        step(k) = -step(k) - step_sign;
      }
    }
  }
  // order the arrays according to sq_norm
  std::vector<Eigen::Index> i_sort = Argsort(sq_norm);

  // TODO: ratio test
  a = Z.inverse() * z.block(0, i_sort.at(0), n, 1);

#ifdef LAMBDA_DEBUG
  std::cout << "z" << z << std::endl;
  std::cout << "square_norm" << sq_norm << std::endl;
  // std::cout << "a" << a << std::endl;
#endif // LAMBDA_DEBUG

}

Eigen::Index PBD_Lambda::Argmax(const Eigen::VectorXd& x)
{
  Eigen::Index row;
  x.maxCoeff(&row);
  return row;
}

Eigen::Index PBD_Lambda::Argmin(const Eigen::VectorXd& x)
{
  Eigen::Index row;
  x.minCoeff(&row);
  return row;
}

std::vector<Eigen::Index> PBD_Lambda::Argsort(Eigen::VectorXd x)
{
  Eigen::VectorXd x_cpy = x;
  size_t n = x.size();
  Eigen::Index i_min;
  std::vector<Eigen::Index> indeces;
  while (indeces.size() < n)
  {
    i_min = Argmin(x_cpy);
    indeces.push_back(i_min);
    RemoveRow(x_cpy, i_min);
  }
  return indeces;
}

void PBD_Lambda::RemoveRow(Eigen::VectorXd& vector, unsigned int row)
{
  // Á‚µ‚Ä‚¢‚­‚Æ‚¸‚ê‚é‚©‚çŒã‚ë‚©‚ç
  // for (int row = end_row; row >= begin_row; --row) {
    unsigned int numRows = vector.rows() - 1;

    if (row < numRows)
      vector.block(row, 0, numRows - row, 1) = vector.bottomRows(numRows - row);

    vector.conservativeResize(numRows, 1);
  // }
}


