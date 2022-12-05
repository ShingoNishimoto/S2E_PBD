#include "PBD_const.h"
#include "PBD_Lambda.h"
#include <iostream>
#include <cassert>
#define _USE_MATH_DEFINES
#include <math.h>
#include "../../Library/VectorTool.hpp"

#define NUM_SOL_CANDIDATE (2)
static const double ratio_threshold = 3.0;

#define LAMBDA_DEBUG

PBD_Lambda::PBD_Lambda(vector<Ambiguity*> N_est_vec, vector<Eigen::MatrixXd>& P_N, const vector<vector<int>> observed_gnss_sat_ids): observed_gnss_sat_ids_(observed_gnss_sat_ids), vec_N_(N_est_vec)
{
  P_ddN_ = InitializeCovariance(P_N);
  n_ = N_hat_.rows();
  Z_ = Eigen::MatrixXd::Identity(n_, n_);
  Eigen::LDLT<Eigen::MatrixXd> LDLTOfP_ddN(P_ddN_.transpose());
  L_ = LDLTOfP_ddN.matrixL();
  D_ = LDLTOfP_ddN.vectorD(); // ここで負の値が出るのは数値誤差的にありえる．
  z_ = Eigen::MatrixXd::Zero(n_, NUM_SOL_CANDIDATE);
  sq_norm_ = Eigen::VectorXd::Zero(NUM_SOL_CANDIDATE);

#ifdef LAMBDA_DEBUG
  std::cout << "P_ddN" << P_ddN_ << std::endl;
  std::cout << "N_est" << N_hat_ << std::endl;
#endif // LAMBDA_DEBUG
  // N_ ?
  // L_.triangularView<Eigen::UnitLower>() = Eigen::MatrixXd::Identity(n, n);
}

PBD_Lambda::~PBD_Lambda() {}

bool PBD_Lambda::Solve(void)
{
  // 部分的に解けている場合にどうするべきなのかわからん．
  if (std::find(vec_N_.at(0)->is_fixed.begin(), vec_N_.at(0)->is_fixed.end(), false) == vec_N_.at(0)->is_fixed.end() &&
      std::find(vec_N_.at(1)->is_fixed.begin(), vec_N_.at(1)->is_fixed.end(), false) == vec_N_.at(1)->is_fixed.end()) return true; // 全部解けているときは何もしない．

  Reduction();
  Eigen::VectorXd z_hat = Z_.transpose() * N_hat_;
  Eigen::MatrixXd P_z_hat = Z_.transpose() * P_ddN_ * Z_; // PARじゃないとつかわない．
  vector<Eigen::VectorXd> N_candidates = IntegerLeastSquares(NUM_SOL_CANDIDATE, z_hat);

  // if (!SuccessRateTest(N_candidates)) return false;
  if (!DiscriminationTest(N_candidates)) return false;
#ifdef LAMBDA_DEBUG
  std::cout << "D_" << D_ << std::endl;
  std::cout << "Z_" << Z_ << std::endl;
  std::cout << "L_" << L_ << std::endl;
  std::cout << "N_hat_" << N_hat_ << std::endl;
  std::cout << "N_" << N_ << std::endl;
#endif // LAMBDA_DEBUG
  RecoverNoDifference();

  // シンプルにここでtestを入れてない影響で，違った解に収束しそれによって新しい衛星に対して収束していない可能性がある.
  return true;
}

// 参照衛星の選択 -> 二重差分の計算部分を実装する．
Eigen::MatrixXd PBD_Lambda::InitializeCovariance(vector<Eigen::MatrixXd> P_N)
{
  // まず，ΔNに対応するPを計算する．その中で最も誤差が小さいものを参照衛星とする（ただ今の状況だとどれもほぼ同じ値になってそう．．）．
  const int n_common = observed_gnss_sat_ids_.at(2).size();
  Eigen::VectorXd d_N_est = Eigen::VectorXd::Zero(n_common);
  Eigen::MatrixXd M_sd = Eigen::MatrixXd::Zero(n_common, 2*n_common);
  Eigen::MatrixXd P_N_all = Eigen::MatrixXd::Zero(2*n_common, 2*n_common); // 一旦受信機間の相関はないと考えてよさそう．

  // ここ，サイズずれるとだるい．
  P_N_all.topLeftCorner(n_common, n_common) = P_N.at(0).topLeftCorner(n_common, n_common);
  P_N_all.bottomRightCorner(n_common, n_common) = P_N.at(1).topLeftCorner(n_common, n_common);
  // commonから探索する．
  for (int ch = 0; ch < observed_gnss_sat_ids_.at(2).size(); ch++)
  {
    const int gnss_id = observed_gnss_sat_ids_.at(2).at(ch);
    // SDCP計算するときにまとめてやりたい気持ちはあるな，chを記録しておけばいいだけな気はする．
    const int main_ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(0), gnss_id);
    const int target_ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(1), gnss_id);
    d_N_est(ch) = vec_N_.at(1)->N.at(target_ch) - vec_N_.at(0)->N.at(main_ch);
    M_sd(main_ch, main_ch) = -1.0; M_sd(target_ch, n_common + target_ch) = 1.0;
    if (main_ch >= n_common || target_ch >= n_common)
    {
      // Pの方を入れ替える操作をしないといけない．一旦なさそうなので無視する？
      abort();
    }
  }
  Eigen::MatrixXd P_dN = M_sd * P_N_all * M_sd.transpose();

  // まず最も分散の小さな衛星を選択する．
  double min_var = 1e5; // 大きめの値．
  int min_index;
  for (int i = 0; i < P_dN.rows(); i++)
  {
    if (P_dN(i, i) < min_var)
    {
      min_var = P_dN(i, i); min_index =i;
    }
  }
  master_d_N_ = std::make_pair(min_index, std::round(d_N_est(min_index))); // 最初は丸める．

  // その次に二重差分Nと誤差共分散を計算．
  Eigen::MatrixXd M_dd = Eigen::MatrixXd::Zero(n_common - 1, n_common);
  Eigen::MatrixXd P_ddN = Eigen::MatrixXd::Zero(n_common - 1, n_common - 1);
  int pos = 0;
  for (int i = 0; i < d_N_est.rows(); i++)
  {
    if (i == min_index) continue;
    M_dd(pos, pos) = 1.0; M_dd(pos, min_index) = -1.0;
    pos++;
  }
  N_hat_ = M_dd * d_N_est;
  P_ddN = M_dd * P_dN * M_dd.transpose();
  return P_ddN;
}


void PBD_Lambda::IntegerGaussTransformation(const int i, const int j)
{
  int mu;
  mu = (int)round(L_(i, j));
  if (mu != 0)
  {
    L_.block(i, j, n_ - i, 1) = L_.block(i, j, n_ - i, 1) - mu * L_.block(i, i, n_ - i, 1);
    Z_.block(0, j, n_, 1) = Z_.block(0, j, n_, 1) - mu * Z_.block(0, i, n_, 1);
    N_hat_(j) = N_hat_(j) - mu * N_hat_(i);
  }

#ifdef LAMBDA_DEBUG
  // std::cout << "L_" << L_ << std::endl;
  // std::cout << "Z_" << Z_ << std::endl;
  // std::cout << "N_" << N_hat_ << std::endl;
  // std::cout << "LZ" << L_*Z_ << std::endl;
#endif // LAMBDA_DEBUG
}

// kの許される範囲 1 < k < n-3
void PBD_Lambda::Permute(const int k, const double delta)
{
  double ita = D_(k) / delta;
  double lambda = D_(k + 1) * L_(k + 1, k) / delta;
  D_(k) = ita * D_(k + 1);
  D_(k + 1) = delta;

  Eigen::Matrix2d temp;
  temp << -L_(k + 1, k), 1,
            ita, lambda;
  if (k > 0) L_.block(k, 0, 2, k) = temp * L_.block(k, 0, 2, k);
  L_(k + 1, k) = lambda;

  // swap
  Eigen::MatrixXd L_temp = L_.block(k + 2, k, n_ - (k + 2), 1);
  L_.block(k + 2, k, n_ - (k + 2), 1) = L_.block(k + 2, k + 1, n_ - (k + 2), 1);
  L_.block(k + 2, k + 1, n_ - (k + 2), 1) = L_temp;

  Z_.col(k).swap(Z_.col(k + 1));
  // N_hat_.row(k).swap(N_hat_.row(k+1));
  // // Pのswapは言及されてないが，ratioテストで必要なので．
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

  // 前処理．本当はなくていいはずだが，数値誤差を考えるとあり得るので．
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
        IntegerGaussTransformation(i, k); // Lの相関部分が0になってしまったら何も進まなくなる．
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

vector<Eigen::VectorXd> PBD_Lambda::IntegerLeastSquares(const int n_cands, Eigen::VectorXd z_hat)
{
  assert(z_hat.size() > 0);

  double chi2 = 1.0 * 10e+18;
  Eigen::VectorXd dist = Eigen::VectorXd::Zero(n_);
  int search_finished = 0;
  int count = -1; // the number of candidates

  Eigen::VectorXd z_cond = Eigen::VectorXd::Zero(n_);
  Eigen::VectorXd z_round = Eigen::VectorXd::Zero(n_);
  Eigen::VectorXi step = Eigen::VectorXi::Zero(n_);

  z_cond(n_ - 1) = z_hat(n_ - 1);
  z_round(n_ - 1) = round(z_cond(n_ - 1));

  double diff = z_cond(n_ - 1) - z_round(n_ - 1);
  step(n_ - 1) = (diff >= 0) ? 1 : -1;

  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n_, n_);
  int k = n_ - 1;
  int i_max = n_cands - 1; // Initially, maximum F(z_) is at n_cands

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
        if (count < n_cands - 2)
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
      if (k == n_ - 1) search_finished = 1;
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
  // order the arrays according to sq_norm_
  vector<Eigen::Index> i_sort = Argsort(sq_norm_);

  vector<Eigen::VectorXd> N_candidates_;
  for (int j = 0; j < NUM_SOL_CANDIDATE; j++)
  {
    // 初めにZ_.transpose()をかけたので普通のZ_をかける．
    N_candidates_.push_back(Z_ * z_.block(0, i_sort.at(j), n_, 1));
  }

#ifdef LAMBDA_DEBUG
  std::cout << "z_" << z_ << std::endl;
  std::cout << "square_norm" << sq_norm_ << std::endl;
  std::cout << "N_1" << N_candidates_.at(0) << std::endl;
  std::cout << "N_2" << N_candidates_.at(1) << std::endl;
#endif // LAMBDA_DEBUG
  return N_candidates_; // 確度の高い最初の二つだけは解けてそう．そのほかがバグっている？そもそもこの精度では解けない？
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
  // 参照渡しじゃないからこれ要らん？
  Eigen::VectorXd x_cpy = x;
  size_t n = x.size();
  Eigen::Index i_min;
  std::vector<Eigen::Index> indexes;
  while (indexes.size() < n)
  {
    i_min = Argmin(x_cpy);
    indexes.push_back(i_min);
    RemoveRow(x_cpy, i_min);
  }
  return indexes;
}

void PBD_Lambda::RemoveRow(Eigen::VectorXd& vector, unsigned int row)
{
  // 消していくとずれるから後ろから
  // for (int row = end_row; row >= begin_row; --row) {
    unsigned int numRows = vector.rows() - 1;

    if (row < numRows)
      vector.block(row, 0, numRows - row, 1) = vector.bottomRows(numRows - row);

    vector.conservativeResize(numRows, 1);
  // }
}

bool PBD_Lambda::SuccessRateTest(vector<Eigen::VectorXd> N_candidates)
{
  // todo
  return true;
}

// ratio testともいう．
bool PBD_Lambda::DiscriminationTest(vector<Eigen::VectorXd> N_candidates)
{
  // 1番目と2番目の解を比較する．
  vector<double> ratio;
  for (int i = 0; i < NUM_SOL_CANDIDATE; i++)
  {
    ratio.push_back((N_hat_ - N_candidates.at(i)).transpose() * P_ddN_.inverse() * (N_hat_ - N_candidates.at(i)));
  }

  if ((ratio.at(1) / ratio.at(0)) > ratio_threshold) return true;
  return false;
}

void PBD_Lambda::RecoverNoDifference(void)
{
  vector<int> sd_N{master_d_N_.second};
  for (int i = 0; i < n_; i++)
  {
    sd_N.push_back(N_(i));
  }

  // originalに置いて精度はまあまあでいいから丸めて差分の不定性が正しければいい？<-真値とはほぼ遠いがいいのか？
  for (int i = 0; i < observed_gnss_sat_ids_.at(0).size(); i++)
  {
    vec_N_.at(0)->N.at(i) = std::round(vec_N_.at(0)->N.at(i));
    vec_N_.at(0)->is_fixed.at(i) = true;
  }

  // この部分は関数かすべきか？
  for (int ch = 0; ch < observed_gnss_sat_ids_.at(2).size(); ch++)
  {
    const int gnss_id = observed_gnss_sat_ids_.at(2).at(ch);
    const int main_ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(0), gnss_id);
    const int target_ch = GetIndexOfStdVector(observed_gnss_sat_ids_.at(1), gnss_id);
    vec_N_.at(1)->N.at(target_ch) = sd_N.at(ch) + vec_N_.at(0)->N.at(main_ch); // 丸めた値で足す．
    vec_N_.at(1)->is_fixed.at(target_ch) = true;
  }
}

