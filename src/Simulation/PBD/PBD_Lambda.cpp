#include "PBD_const.h"
#include "PBD_Lambda.h"
#include <iostream>
#include <cassert>
#define _USE_MATH_DEFINES
#include <math.h>
#include "../../Library/VectorTool.hpp"

#define NUM_SOL_CANDIDATE (2)
static const double ratio_threshold = 3.0;

// #define LAMBDA_DEBUG

PBD_Lambda::PBD_Lambda(vector<Ambiguity*> N_est_vec, vector<Eigen::MatrixXd>& P_N, const vector<vector<int>> observed_gnss_sat_ids, const int num_single_state): observed_gnss_sat_ids_(observed_gnss_sat_ids), vec_N_(N_est_vec)
{
  P_ddN_ = InitializeCovariance(P_N, num_single_state);

// debug ++++++++++++++++++
  // P_ddN_ = Eigen::Matrix3d::Zero();
  // P_ddN_ << 6.2900, 5.9780, 0.5440,
  //           5.9780, 6.2920, 2.3400,
  //           0.5440, 2.3400, 6.2880;
  // N_hat_ = Eigen::Vector3d::Zero();
  // N_hat_ << 5.4500, 3.1000, 2.9700;
// ++++++++++++++++++++++++

  n_ = N_hat_.rows();
  Z_ = Eigen::MatrixXd::Identity(n_, n_);
  L_ = Eigen::MatrixXd::Zero(n_, n_);
  D_ = Eigen::VectorXd::Zero(n_);
  z_ = Eigen::MatrixXd::Zero(n_, NUM_SOL_CANDIDATE);
  sq_norm_ = Eigen::VectorXd::Zero(NUM_SOL_CANDIDATE);
  N_ = N_hat_;

#ifdef LAMBDA_DEBUG
  std::cout << "P_ddN" << P_ddN_ << std::endl;
  std::cout << "D" << D_ << std::endl;
  std::cout << "L" << L_ << std::endl;
  std::cout << "N_est" << N_hat_ << std::endl;
#endif // LAMBDA_DEBUG
  // N_ ?
  // L_.triangularView<Eigen::UnitLower>() = Eigen::MatrixXd::Identity(n, n);
}

PBD_Lambda::~PBD_Lambda() {}

bool PBD_Lambda::Solve(METHOD method)
{
  // �����I�ɉ����Ă���ꍇ�ɂǂ�����ׂ��Ȃ̂��킩���D
  if (std::find(vec_N_.at(0)->is_fixed.begin(), vec_N_.at(0)->is_fixed.end(), false) == vec_N_.at(0)->is_fixed.end() &&
      std::find(vec_N_.at(1)->is_fixed.begin(), vec_N_.at(1)->is_fixed.end(), false) == vec_N_.at(1)->is_fixed.end()) return true; // �S�������Ă���Ƃ��͉������Ȃ��D

  // ����l�s��ɂȂ��Ă��邩�m�F�D
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES(P_ddN_);
  Eigen::VectorXd e = ES.eigenvalues();
  if ((e.array() < 0).any()) return false; // �����I�ɉ������s��͂����Ŗ���͂�����Ă��܂��D

  // Eigen::LDLT<Eigen::MatrixXd> LDLTOfP_ddN(P_ddN_);
  // L_ = LDLTOfP_ddN.matrixL();
  // D_ = LDLTOfP_ddN.vectorD(); // �����ŕ��̒l���o��̂͐��l�덷�I�ɂ��肦��D
  RobustCholeskyDecomposition(P_ddN_);

  Reduction();

#ifdef LAMBDA_DEBUG
  std::cout << "D_" << D_ << std::endl;
  std::cout << "Z_" << Z_ << std::endl;
  std::cout << "L_" << L_ << std::endl;
  std::cout << "N_hat_" << N_hat_ << std::endl;
#endif // LAMBDA_DEBUG

  // �����Ő������Ə������ɕ�������D
  Eigen::VectorXi int_parts = N_hat_.cast<int>();
  Eigen::VectorXd float_parts = N_hat_ - int_parts.cast<double>();
  Eigen::VectorXd z_hat = Z_.transpose() * float_parts;
  Eigen::MatrixXd P_z_hat = Z_.transpose() * P_ddN_ * Z_; // PAR����Ȃ��Ƃ���Ȃ��D

  // if (!(CalculateSuccessRate(D_) > 0.9)) return false;

  int fixed_num;
  switch (method)
  {
  case METHOD::ILS:
    if (IntegerLeastSquares(z_hat)) fixed_num = n_;
    else fixed_num = 0;
    break;

  case METHOD::PAR_ILS:
    // const double p0 = 0.8; // �����ƌ��߂�D
    fixed_num = PartialSearch(z_hat, P_z_hat, 0.8);
    break;

  default:
    fixed_num = 0;
    break;
  }

  if (!fixed_num)
  {
    RecoverNoDifference(fixed_num);
    return false;
  }

  // order the arrays according to sq_norm_
  vector<Eigen::Index> i_sort = Argsort(sq_norm_);
  vector<Eigen::VectorXd> N_candidates;
  for (int j = 0; j < NUM_SOL_CANDIDATE; j++)
  {
    // ���߂ɐ������������āCZ_.transpose()���������̂ŁC�����߂��D
    N_candidates.push_back((Z_.transpose()).inverse() * z_.block(0, i_sort.at(j), n_, 1) + int_parts.cast<double>());
  }

#ifdef LAMBDA_DEBUG
  std::cout << "z_" << z_ << std::endl;
  std::cout << "square_norm" << sq_norm_ << std::endl;
  std::cout << "N_1" << N_candidates.at(0) << std::endl;
  std::cout << "N_2" << N_candidates.at(1) << std::endl;
#endif // LAMBDA_DEBUG

  // if (!DiscriminationTest(N_candidates))
  // {
  //   vec_N_.at(0)->is_fixed.assign(vec_N_.at(0)->is_fixed.size(), false);
  //   vec_N_.at(1)->is_fixed.assign(vec_N_.at(1)->is_fixed.size(), false);
  //   return false;
  // }

  N_ = N_candidates.at(0);
  RecoverNoDifference(fixed_num);

  return true;
}

// �Q�Ɖq���̑I�� -> ��d�����̌v�Z��������������D
Eigen::MatrixXd PBD_Lambda::InitializeCovariance(vector<Eigen::MatrixXd> P_N, const int num_single_state)
{
  // �ϐ�����
  const int num_double_state = 2 * num_single_state;
  const int n_main = observed_gnss_sat_ids_.at(0).size();
  const int n_target = observed_gnss_sat_ids_.at(1).size();
  const int n_common = observed_gnss_sat_ids_.at(2).size();
  const int all_state_num = num_double_state + n_main + n_target;
  const int sd_state_num = num_double_state + n_common;

  // �܂��C��N�ɑΉ�����P���v�Z����D
  Eigen::VectorXd d_N_est = Eigen::VectorXd::Zero(n_common);
  Eigen::MatrixXd M_sd = Eigen::MatrixXd::Zero(n_common, n_main + n_target);
  Eigen::MatrixXd P_N_all = Eigen::MatrixXd::Zero(n_main + n_target, n_main + n_target); // ��U��M�@�Ԃ̑��ւ͂Ȃ��ƍl���Ă悳�����D
  Eigen::MatrixXd M1 = Eigen::MatrixXd::Zero(sd_state_num, all_state_num); // ���ԍs��
  M1.topLeftCorner(num_single_state, num_single_state) = Eigen::MatrixXd::Identity(num_single_state, num_single_state);
  M1.block(num_single_state, num_single_state + n_main, num_single_state, num_single_state) = Eigen::MatrixXd::Identity(num_single_state, num_single_state);

  // �����C�T�C�Y�����Ƃ��邢�D
  P_N_all.topLeftCorner(n_common, n_common) = P_N.at(0).topLeftCorner(n_common, n_common);
  P_N_all.bottomRightCorner(n_common, n_common) = P_N.at(1).topLeftCorner(n_common, n_common);
  // common����T������D
  for (int ch = 0; ch < n_common; ch++)
  {
    const int gnss_id = observed_gnss_sat_ids_.at(2).at(ch);
    // SDCP�v�Z����Ƃ��ɂ܂Ƃ߂Ă�肽���C�����͂���ȁCch���L�^���Ă����΂��������ȋC�͂���D
    const int main_ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(0), gnss_id);
    const int target_ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(1), gnss_id);
    d_N_est(ch) = vec_N_.at(1)->N.at(target_ch) - vec_N_.at(0)->N.at(main_ch);
    M_sd(ch, main_ch) = -1.0; M_sd(ch, n_common + target_ch) = 1.0;
    if (main_ch >= n_common || target_ch >= n_common)
    {
      // P�̕������ւ��鑀������Ȃ��Ƃ����Ȃ��D��U�Ȃ������Ȃ̂Ŗ�������H
      abort();
    }
  }
  Eigen::MatrixXd P_dN = M_sd * P_N_all * M_sd.transpose();
  M1.block(num_double_state, num_single_state, n_common, n_main) = M_sd.topLeftCorner(n_common, n_main);
  M1.block(num_double_state, num_double_state + n_main, n_common, n_target) = M_sd.bottomRightCorner(n_common, n_target);

  // �܂��ł����U�̏����ȉq�����Q�Ɖq���Ƃ���D
  // 2��ڈȍ~�͉����Ă�������x�[�X�Ɏ��{���Ă����ׂ��D
  int ref_index;
  const int ref_index_candidate = GetIndexOfStdVector(vec_N_.at(0)->is_fixed, true);
  if (ref_index_candidate >= 0 && ref_index_candidate < n_common) // vec_N_.at(0)->is_fixed.size()
  {
    ref_index = ref_index_candidate;
  }
  else
  {
    double min_var = 1e5; // �傫�߂̒l�D
    for (int i = 0; i < P_dN.rows(); i++)
    {
      if (P_dN(i, i) < min_var)
      {
        min_var = P_dN(i, i); ref_index =i;
      }
    }
  }

  master_d_N_ = std::make_pair(ref_index, std::round(d_N_est(ref_index))); // �ŏ��͊ۂ߂�Dpull-in�̈�͂���ł����̂��H

  // ���̎��ɓ�d����N�ƌ덷�����U���v�Z�D
  Eigen::MatrixXd M2 = Eigen::MatrixXd::Zero(sd_state_num - 1, sd_state_num);
  M2.topLeftCorner(num_double_state, num_double_state) = Eigen::MatrixXd::Identity(num_double_state, num_double_state);
  Eigen::MatrixXd M_dd = Eigen::MatrixXd::Zero(n_common - 1, n_common);
  Eigen::MatrixXd P_ddN = Eigen::MatrixXd::Zero(n_common - 1, n_common - 1);
  int pos = 0;
  for (int i = 0; i < d_N_est.rows(); i++)
  {
    if (i == ref_index) continue;

    M_dd(pos, i) = 1.0; M_dd(pos, ref_index) = -1.0;
    gnss_ids_.push_back(observed_gnss_sat_ids_.at(2).at(i)); // �ŏ��̕��т��L�^�D
    pos++;
  }
  N_hat_ = M_dd * d_N_est;
  P_ddN = M_dd * P_dN * M_dd.transpose();

  M2.bottomRightCorner(n_common - 1, n_common) = M_dd;
  M2M1_ = M2 * M1;

#ifdef LAMBDA_DEBUG
  std::cout << "P_N0" << P_N.at(0) << std::endl;
  std::cout << "P_N1" << P_N.at(1) << std::endl;
  std::cout << "P_dN" << P_dN << std::endl;
  std::cout << "P_ddN" << P_ddN << std::endl;
#endif // LAMBDA_DEBUG

  return P_ddN;
}

// LDLt decomposition �����ł�reordering�͓����ĂȂ��D
void PBD_Lambda::RobustCholeskyDecomposition(Eigen::MatrixXd P)
{
  const int n = P.rows();

  for (int i = n -1; i >= 0; i--)
  {
    D_(i) = P(i, i); // �ǂ̍s���J�n�_�ɂ��邩���C��t����Ζ��͋N����Ȃ��͂������C������������Ȃ��Ƃ��l���邱�Ǝ��̊Ԉ���Ă���?
    // if (D_(i) <= 0)
    L_.block(i, 0, 1, i + 1) = P.block(i, 0, 1, i + 1).array() / sqrt(P(i, i));

    for (int j = 0; j < i; j++) P.block(j, 0, 1, j + 1) -= L_.block(i, 0, 1, j + 1) * L_(i, j);

    L_.block(i, 0, 1, i + 1) /= L_(i, i);
  }

  if ((D_.array() < 1e-10).any())
  {
    abort();
  }
}

void PBD_Lambda::IntegerGaussTransformation(const int i, const int k)
{
  int mu;
  mu = (int)round(L_(i, k));
  if (mu != 0)
  {
    L_.block(i, k, n_ - i, 1) -= (mu * L_.block(i, i, n_ - i, 1).array()).matrix();
    Z_.col(i) += mu * Z_.col(k);
    // N_hat_(k) = N_hat_(k) - mu * N_hat_(i);
  }

#ifdef LAMBDA_DEBUG
  // std::cout << "L_" << L_ << std::endl;
  // std::cout << "Z_" << Z_ << std::endl;
  // std::cout << "N_" << N_hat_ << std::endl;
  // std::cout << "LZ" << L_*Z_ << std::endl;
#endif // LAMBDA_DEBUG
}

// k�̋������͈� 1 < k < n-3
void PBD_Lambda::Permute(const int k, const double delta)
{
  double ita = D_(k) / delta;
  double lambda = D_(k + 1) * L_(k + 1, k) / delta;
  D_(k) = ita * D_(k + 1);
  D_(k + 1) = delta;

  Eigen::Matrix2d temp;
  temp << -L_(k + 1, k), 1,
            ita, lambda;
  // if (k > 0)
  L_.block(k, 0, 2, k) = temp * L_.block(k, 0, 2, k);
  L_(k + 1, k) = lambda;

  // swap
  Eigen::MatrixXd L_temp = L_.block(k + 2, k, n_ - (k + 2), 1);
  L_.block(k + 2, k, n_ - (k + 2), 1) = L_.block(k + 2, k + 1, n_ - (k + 2), 1);
  L_.block(k + 2, k + 1, n_ - (k + 2), 1) = L_temp;

  Z_.col(k).swap(Z_.col(k + 1));
  // ���̓���ւ��ɂ��gnss_id�̑Ή����L�^
  std::swap(gnss_ids_.at(k), gnss_ids_.at(k+1));
  // N_hat_.row(k).swap(N_hat_.row(k+1));
  // // P��swap�͌��y����ĂȂ����Cratio�e�X�g�ŕK�v�Ȃ̂ŁD
  // P_ddN_.row(k).swap(P_ddN_.row(k+1));
  // P_ddN_.col(k).swap(P_ddN_.col(k+1));

#ifdef LAMBDA_DEBUG
  // std::cout << "D_" << D_ << std::endl;
  // std::cout << "Z_" << Z_ << std::endl;
  // std::cout << "N_" << N_hat_ << std::endl;
#endif // LAMBDA_DEBUG
}

void PBD_Lambda::Reduction(void)
{
  int k, l;
  k = n_ - 2;
  l = n_ - 2;
  Eigen::VectorXd D_bar = Eigen::VectorXd::Zero(n_);

  // �O�����D�{���͂Ȃ��Ă����͂������C���l�덷���l����Ƃ��蓾��̂ŁD
  for (int i = 0; i < D_.rows(); i++)
  {
    if (D_(i) < 0) D_(i) *= -1;
  }

  while (k >= 0)
  {
    assert(k < n_ - 1);
    if (k <= l)
    {
      for (int i = k + 1; i < n_; i++)
      {
        IntegerGaussTransformation(i, k); // L�̑��֕�����0�ɂȂ��Ă��܂����牽���i�܂Ȃ��Ȃ�D
      }
    }
    D_bar(k + 1) = D_(k) + pow(L_(k + 1, k), 2) * D_(k + 1);
    if (D_bar(k + 1) < D_(k + 1))
    {
      Permute(k, D_bar(k + 1));
      l = k;
      k = n_ - 2;
    }
    else
    {
      k--;
    }
  }

  // calculate Z transformation matrix
  Eigen::MatrixXd Z_float = (Z_.transpose()).inverse();
  for (int i = 0; i < n_; i++)
  {
    for (int j = 0; j < n_; j++) Z_(i, j) = round(Z_float(i, j));
  }
}

bool PBD_Lambda::IntegerLeastSquares(Eigen::VectorXd z_hat) // const int c_cands
{
  const int n = z_hat.rows();
  assert(n > 0);

  double chi2 = 1.0 * 10e+18;
  Eigen::VectorXd dist = Eigen::VectorXd::Zero(n);
  int search_finished = 0;
  int count = -1; // the number of candidates

  z_ = Eigen::MatrixXd::Zero(n, NUM_SOL_CANDIDATE);
  Eigen::VectorXd z_cond = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd z_round = Eigen::VectorXd::Zero(n);
  Eigen::VectorXi step = Eigen::VectorXi::Zero(n);

  z_cond(n - 1) = z_hat(n - 1);
  z_round(n - 1) = round(z_cond(n - 1));

  double diff = z_cond(n - 1) - z_round(n - 1);
  step(n - 1) = (diff >= 0) ? 1 : -1;

  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
  int k = n - 1;
  int i_max = NUM_SOL_CANDIDATE - 1; // Initially, maximum F(z_) is at n_cands

  while (!search_finished)
  {
    double dist_new = dist(k) + pow(diff, 2) / D_(k);
    if (dist_new < chi2)
    {
      if (k != 0)
      {
        k--;
        dist(k) = dist_new;
        S.block(k, 0, 1, k) = S.block(k + 1, 0, 1, k) + (z_round(k + 1) - z_cond(k + 1)) * L_.block(k + 1, 0, 1, k);

        z_cond(k) = z_hat(k) + S(k, k);
        z_round(k) = round(z_cond(k));
        diff = z_cond(k) - z_round(k);
        step(k) = (diff >= 0) ? 1 : -1;
      }
      else
      {
        // store the found candidate and try next valid integer
        if (count < NUM_SOL_CANDIDATE - 2)
        {
          count++;
          z_.col(count) = z_round;
          sq_norm_(count) = dist_new;
        }
        else
        {
          z_.col(i_max) = z_round;
          sq_norm_(i_max) = dist_new;

          i_max = Argmax(sq_norm_);
          chi2 = sq_norm_(i_max);
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

  // ������ratio test�����΂��������������D����Ȃ��͑I�ׂ�������������
  vector<Eigen::Index> i_sort = Argsort(sq_norm_);
  if (sq_norm_(i_sort.at(1)) / sq_norm_(i_sort.at(0)) < ratio_threshold)
  {
    z_ = z_hat; // �����̂܂�
    return false;
  }

  return true;
}

int PBD_Lambda::PartialSearch(Eigen::VectorXd z_hat, Eigen::MatrixXd P_z, const double p0)
{
  double p_s = CalculateSuccessRate(D_);

  int k = 0;
  while (p_s < p0 && k < n_ - 1)
  {
    k++;
    // compute the bootstrapped success rate if the last n-k+1 ambiguities would be fixed
    p_s = CalculateSuccessRate(D_.bottomRows(n_ - k));
  }

  if (p_s < p0) return 0; // ������������݂��Ȃ��D�Ԃ�l���l����D

  if (!IntegerLeastSquares(z_hat.bottomRows(n_ - k)))
  {
    z_ = z_hat;
    return 0;
  }

  Eigen::MatrixXd z_par = z_;
  Eigen::MatrixXd P_z_par = P_z.bottomRightCorner(n_ - k, n_ - k);
  // Eigen::MatrixXd Z_par = Z_.block(0, k, n_, n_ - k);

  // first k-1 ambiguities are adjusted based on correlation with the fixed ambiguities
  Eigen::MatrixXd QP = P_z.block(0, k, k, n_ - k) * P_z_par.inverse();
  Eigen::MatrixXd z_float = Eigen::MatrixXd::Zero(k, NUM_SOL_CANDIDATE);
  for (int i = 0; i < NUM_SOL_CANDIDATE; i++)
  {
    z_float.col(i) = z_hat.topRows(k) - QP * (z_hat.bottomRows(n_ - k) - z_par.col(i));
  }

  z_ = Eigen::MatrixXd::Zero(n_, NUM_SOL_CANDIDATE); // ������
  z_.topRows(k) = z_float;
  z_.bottomRows(n_ - k) = z_par;

#ifdef LAMBDA_DEBUG
  std::cout << "z_" << z_ << std::endl;
#endif // LAMBDA_DEBUG

  return n_ - k; // fixed_num
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
  size_t n = x.size();
  Eigen::Index i_min;
  std::vector<Eigen::Index> indexes;
  while (indexes.size() < n)
  {
    i_min = Argmin(x);
    indexes.push_back(i_min);
    // RemoveRow(x, i_min);
    // �ǂ��Ȃ����ǉ��}���u�Ƃ��Č��̏ꏊ�ɑ傫�Ȓl������D
    x(i_min) = 1e18;
  }
  return indexes;
}

void PBD_Lambda::RemoveRow(Eigen::VectorXd& vector, unsigned int row)
{
  // �����Ă����Ƃ���邩���납��
  // for (int row = end_row; row >= begin_row; --row) {
    unsigned int numRows = vector.rows() - 1;

    if (row < numRows)
      vector.block(row, 0, numRows - row, 1) = vector.bottomRows(numRows - row);

    vector.conservativeResize(numRows, 1);
  // }
}

double PBD_Lambda::CalculateSuccessRate(Eigen::VectorXd D)
{
  double prob = 1.0;
  const int n = D.rows();
  for (int i = 0; i < n; i++)
  {
    // prob *= pow(1 - exp(-0.5 * 0.25/D(i)), 0.5);
    prob *= erf(0.5/sqrt(D(i)) / sqrt(2.0)); // �ǂ����ł��ς���ȁD
  }

  return prob;
}

// ratio test�Ƃ������D
bool PBD_Lambda::DiscriminationTest(vector<Eigen::VectorXd> N_candidates) // �����fix�������̂ɑ΂��Ă������{���ׂ��Ȃ̂ł́H
{
  // 1�Ԗڂ�2�Ԗڂ̉����r����D
  vector<double> ratio;
  for (int i = 0; i < NUM_SOL_CANDIDATE; i++)
  {
    ratio.push_back((N_hat_ - N_candidates.at(i)).transpose() * P_ddN_.inverse() * (N_hat_ - N_candidates.at(i)));
  }

  if ((ratio.at(1) / ratio.at(0)) > ratio_threshold) return true;
  return false;
}

void PBD_Lambda::RecoverNoDifference(const int fixed_num)
{
  vector<double> sd_N;
  int pos = 0;
  for (int i = 0; i < n_ + 1; i++)
  {
    if (i == master_d_N_.first) sd_N.push_back(master_d_N_.second);
    else
    {
      sd_N.push_back(N_(pos) + master_d_N_.second); // N_��double difference�Ȃ̂�single difference�ɖ߂��D
      pos++;
    }
  }

  // gnss_ids_�ɂ�fix���̃\�[�g�̏��ŕ���ł���D�܂艺�̍s�̂��̂���fix���Ă���Ƃ������ƁD
  gnss_ids_.push_back(observed_gnss_sat_ids_.at(2).at(master_d_N_.first));
  std::reverse(gnss_ids_.begin(), gnss_ids_.end());
  for (int i = 0; i < gnss_ids_.size(); i++)
  {
    const int gnss_id = gnss_ids_.at(i);
    const int ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(2), gnss_id);
    const int main_ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(0), gnss_id);
    const int target_ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(1), gnss_id);
    vec_N_.at(1)->N.at(target_ch) = sd_N.at(ch) + vec_N_.at(0)->N.at(main_ch);
    if (i <= fixed_num) // fixed_num��DD�̐��Ȃ̂�N�ōl�����+1
    {
      vec_N_.at(0)->is_fixed.at(main_ch) = true;
      vec_N_.at(1)->is_fixed.at(target_ch) = true;
    }
    else
    {
      vec_N_.at(0)->is_fixed.at(main_ch) = false;
      vec_N_.at(1)->is_fixed.at(target_ch) = false;
    }
  }

#ifdef LAMBDA_DEBUG
  // gnss_ids�̐��ȏ�̕����͂������������ĂȂ��̂�false�ɗ��Ƃ��i�O�̂��߁j�D
  for (int j = gnss_ids_.size(); j < vec_N_.at(0)->N.size(); j++)
  {
    vec_N_.at(0)->is_fixed.at(j) = false;
    vec_N_.at(1)->is_fixed.at(j) = false;
    std::cout << vec_N_.at(0)->N.at(j) << std::endl;
    std::cout << vec_N_.at(1)->N.at(j) << std::endl;
  }
#endif // LAMBDA_DEBUG
}

