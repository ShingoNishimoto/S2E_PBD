#include "PCVEstimation.hpp"
#include <Constant.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Interface/InitInput/IniAccess.h"

#define SH_WLS_DATA_MULTIPLE (10)

PCVEstimation::PCVEstimation(const std::string fname)
{
  IniAccess pcv_conf(fname);
  std::string method_str = pcv_conf.ReadString("PCV", "pcv_estimation_method");
  if (method_str == "SPHERE")
  {
    method_ = PCV_METHOD::SPHERE;
    SphericalHarmonicsInitialization(fname);
  }
  else if (method_str == "ZERNIKE")
  {
    method_ = PCV_METHOD::ZERNIKE;
    // TODO
  }
  else if (method_str == "RESIDUAL")
  {
    method_ = PCV_METHOD::RESIDUAL;
    // TODO
  }
  else
  {
    std::cout << "ERROR: PCV estimation method error!" << std::endl;
  }
}

PCVEstimation::PCVEstimation(){}

PCVEstimation::~PCVEstimation(){}

void PCVEstimation::SetHRaw(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation)
{
  const double azi_rad_i = gnss_observation.GetGnssAzimuthDeg(i) * libra::deg_to_rad;
  const double ele_rad_i = gnss_observation.GetGnssElevationDeg(i) * libra::deg_to_rad;
  const double azi_rad_j = gnss_observation.GetGnssAzimuthDeg(ref_j) * libra::deg_to_rad;
  const double ele_rad_j = gnss_observation.GetGnssElevationDeg(ref_j) * libra::deg_to_rad;
  const Eigen::MatrixXd Hi = CalcLegendreCoeff(azi_rad_i, ele_rad_i) - CalcLegendreCoeff(azi_rad_j, ele_rad_j);
  const int pos = local_pos + V_.rows();
  H_.block(pos, 0, 1, Hi.cols()) = Hi;
}

void PCVEstimation::UpdateReferenceSat(const int count, int& ref_gnss_ch, const double r_sdcp)
{
  // 観測誤差が小さな衛星を参照衛星とする．
  if (r_sdcp < min_variance_)
  {
    min_variance_ = r_sdcp;
    ref_gnss_ch = count;
  }
}

const bool PCVEstimation::Update(const Eigen::VectorXd& V, const Eigen::MatrixXd& W, const double azi_increment, const double ele_increment)
{
  switch (method_)
  {
  case PCV_METHOD::SPHERE:
    return WeightedLeastSquare(V, W, azi_increment, ele_increment);

  case PCV_METHOD::ZERNIKE:
    // TODO
    return false;

  case PCV_METHOD::RESIDUAL:
    // TODO
    return false;

  default:
    return false;
  }
}

void PCVEstimation::SphericalHarmonicsInitialization(const std::string fname)
{
  IniAccess pcv_conf(fname);
  degree_ = pcv_conf.ReadInt("SPHERE", "degree");
  wsl_data_num_ = SH_WLS_DATA_MULTIPLE * (degree_ + 2) * (degree_ + 1);
  // coefficients
  c_.assign(degree_ + 1, vector<double>(degree_ + 1, 0.0));
  s_.assign(degree_ + 1, vector<double>(degree_ + 1, 0.0));
  CS_vec_ = Eigen::VectorXd::Zero((degree_ + 2) * (degree_ + 1));
  InitializeVHW();
}

const Eigen::MatrixXd PCVEstimation::CalcLegendreCoeff(const double azi_rad, const double ele_rad)
{
  x_ = cos(ele_rad) * cos(azi_rad);
  y_ = cos(ele_rad) * sin(azi_rad);
  z_ = sin(ele_rad);

  // Calc V and W
  // int degree_vw = degree_ + 1;
  vector<vector<double>> v(degree_ + 1, vector<double>(degree_ + 1, 0.0));
  vector<vector<double>> w(degree_ + 1, vector<double>(degree_ + 1, 0.0));
  // n_=m_=0
  v[0][0] = 1.0;
  w[0][0] = 0.0;
  m_ = 0;

  while (m_ < degree_) {
    for (n_ = m_ + 1; n_ <= degree_; n_++)
    {
      if (n_ <= m_ + 1)
        v_w_nm_update(&v[n_][m_], &w[n_][m_], v[n_ - 1][m_], w[n_ - 1][m_], 0.0, 0.0);
      else
        v_w_nm_update(&v[n_][m_], &w[n_][m_], v[n_ - 1][m_], w[n_ - 1][m_], v[n_ - 2][m_], w[n_ - 2][m_]);

      // VW(0, m_ * 2 * (degree_ - m_ + 2) + n_)     = v[n_][m_];
      // VW(0, m_ * 2 * (degree_ - m_ + 2) + n_ + 1) = w[n_][m_];
    }
    // next step
    m_++;
    n_ = m_;
    v_w_nn_update(&v[n_][m_], &w[n_][m_], v[n_ - 1][m_ - 1], w[n_ - 1][m_ - 1]);
  }

  // 上の計算時に代入するのは複雑でバグりそうなので避けて別で代入．
  Eigen::MatrixXd VW = Eigen::MatrixXd::Zero(1, (degree_ + 2) * (degree_ + 1));
  VW(0, 0) = v[0][0]; VW(0, 1) = w[0][0];
  for (n_ = 1; n_ <= degree_; n_++)
  {
    for (m_ = 0; m_ <= n_; m_++)
    {
      VW(0, n_ * (n_ + 1) + 2 * m_) = v[n_][m_];
      VW(0, n_ * (n_ + 1) + 2 * m_ + 1) = w[n_][m_];
    }
  }

  return VW;
}


void PCVEstimation::v_w_nn_update(double *v_nn, double *w_nn, const double v_prev, const double w_prev) {
  if (n_ != m_) return;

  double m_d = (double)m_;

  double c;
  c = 2.0 * m_d - 1.0;

  *v_nn = c * (x_ * v_prev - y_ * w_prev);
  *w_nn = c * (x_ * w_prev + y_ * v_prev);
  return;
}

void PCVEstimation::v_w_nm_update(double *v_nm, double *w_nm, const double v_prev, const double w_prev, const double v_prev2, const double w_prev2) {
  if (n_ == m_) return;

  double m_d = (double)m_;
  double n_d = (double)n_;

  double c1 = (2.0 * n_d - 1.0) / (n_d - m_d);
  double c2 = (n_d + m_d - 1.0) / (n_d - m_d);

  *v_nm = c1 * z_ * v_prev - c2 * v_prev2;
  *w_nm = c1 * z_ * w_prev - c2 * w_prev2;
  return;
}

void PCVEstimation::p_n_0_update(double *p_n0, const double p_prev, const double p_prev2)
{
  double n_d = (double)n_;

  *p_n0 = 1.0 / n_d * ((2 * n_d - 1.0) * p_prev - (n_d - 1.0) * p_prev2);
  return;

}

void PCVEstimation::InitializeVHW(void)
{
  const int coef_size = (degree_ + 2) * (degree_ + 1);
  V_ = Eigen::VectorXd::Zero(1);
  H_ = Eigen::MatrixXd::Zero(1, coef_size);
  W_ = Eigen::MatrixXd::Identity(1, 1);// * 1e18; // constraintなので不確定性は0として．

  vector<double> pn0(degree_ + 1);
  pn0[0] = 1.0; pn0[1] = 1.0;
  for (n_ = 0; n_ <= degree_; n_++)
  {
    if (n_ >= 2) p_n_0_update(&pn0[n_], pn0[n_ - 1], pn0[n_ - 2]);
    H_(0, n_ * (n_ + 1)) = pn0[n_];
  }
  // H_ = CalcLegendreCoeff(90 * libra::deg_to_rad, 90 * libra::deg_to_rad);
  // std::cout << "H" << H_ << std::endl;
}

const bool PCVEstimation::WeightedLeastSquare(const Eigen::VectorXd& V, const Eigen::MatrixXd& W, const double azi_increment, const double ele_increment)
{
  const int N = V.rows();
  const int current_size = V_.rows();
  const int new_size = current_size + N;
  const int coef_num = (degree_ + 2) * (degree_ + 1);
  H_.conservativeResize(new_size, coef_num);
  V_.conservativeResize(new_size);
  Eigen::MatrixXd W_cpy = W_;
  W_ = Eigen::MatrixXd::Zero(new_size, new_size);

  // append new matrix and observations
  // H_.block(current_size, 0, N, 3) = H;
  V_.block(current_size, 0, N, 1) = V;
  W_.topLeftCorner(current_size, current_size) = W_cpy;
  W_.block(current_size, current_size, N, N) = W;

  // std::cout << "H" << H_ << std::endl;
  // std::cout << "V" << V_ << std::endl;
  // std::cout << "W" << W_ << std::endl;

  if (new_size >= wsl_data_num_)
  {
    // Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(W_);
    // auto rank = lu_decomp.rank();
    // if (rank < W_.rows()) abort(); // 正則判定．

    // rank落ちしないように制約を課す．
    // Sn0は消してしまう．
    RemoveZeroCols(H_);
    // std::cout << "H" << H_ << std::endl;
    Eigen::VectorXd CS = (H_.transpose() * W_ * H_).inverse() * (H_.transpose() * W_ * V_);

// debug ++++++++++++++++++++++++++++++++++++++++++++++++
    Eigen::VectorXd dd_pcc_after = H_ * CS;
    // std::cout << "H " << H_ << std::endl;
    // std::cout << "after estimated ddpcc" << dd_pcc_after << std::endl;
    // std::cout << "raw ddpcc" << V_ << std::endl;
    double pre_acc = 0;
    for (int i = 0; i < V_.rows(); i++) pre_acc += pow(V_(i), 2.0);
    pre_acc = sqrt(pre_acc) / V_.rows(); // RMS
    // std::cout << "pre accuracy" << pre_acc << std::endl;

    Eigen::VectorXd post_res = V_ - dd_pcc_after;
    double post_acc = 0;
    for (int i = 0; i < post_res.rows(); i++) post_acc += pow(post_res(i), 2.0);
    post_acc = sqrt(post_acc) / post_res.rows();
    // std::cout << "post accuracy" << post_acc << std::endl;
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // 改善した時だけ更新する．
    if (pre_acc > post_acc)
    {
      CS_vec_ = CS;
      SetPcvVecFromSHModel(azi_increment, ele_increment);
      // 収束判定はどうする？
      // if (pre_acc < ddcp_res_thresh && post_acc < ddcp_res_thresh)
      // {
      //   pcv_fixed_ = true;
      // }
      InitializeVHW();
      return true;
    }

    // リセット
    InitializeVHW();
    // V_.resize(0);
    // H_.resize(0, 0);
    // W_.resize(0, 0);
  }
  return false;
}

void PCVEstimation::SetPcvVecFromSHModel(const double azi_increment, const double ele_increment)
{
  // gridに分解してvectorに保存する．
  const int num_azi = (int)(360 / azi_increment);
  const int num_ele = (int)(90 / ele_increment);

  // 格納前にclear
  dpcv_vec_mm_.clear();
  for (int i = 0; i <= num_azi; i++)
  {
    double azimuth = azi_increment * i;
    for (int j = 0; j <= num_ele; j++)
    {
      double elevation = ele_increment * j;
      Eigen::MatrixXd VW = CalcLegendreCoeff(azimuth * libra::deg_to_rad, elevation * libra::deg_to_rad);
      RemoveZeroCols(VW);
      const double dpcv_mm = (VW * CS_vec_)(0);
      dpcv_vec_mm_.push_back(dpcv_mm);
    }
  }
}

// Sn0を消す．
void PCVEstimation::RemoveZeroCols(Eigen::MatrixXd& H)
{
  const int row_num = H.rows();
  for (int i = degree_; i >= 0; i--)
  {
    const int offset = i * (i + 1) + 1;
    const int col_num = H.cols() - 1;
    H.block(0, offset, row_num, col_num - offset) = H.block(0, offset + 1, row_num, col_num - offset);
    H.conservativeResize(row_num, col_num);
  }
}
