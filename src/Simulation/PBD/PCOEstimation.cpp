#include "PCOEstimation.hpp"
#include <Constant.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

#define PCO_WLS_DATA_NUM (1 * 11)

PCOEstimation::PCOEstimation()
{}

PCOEstimation::~PCOEstimation()
{}

void PCOEstimation::SetHVRaw(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp)
{
  const libra::Vector<3> de = gnss_observation.GetGnssDirection_c(i) - gnss_observation.GetGnssDirection_c(ref_j);
  const int pos = local_pos + W_.rows();
  for (int i = 0; i < 3; i++) H_(pos, i) = de[i];
  V_(pos) = res_ddcp;
}

const bool PCOEstimation::CheckDataForEstimation(const int count, int& ref_gnss_ch, const double elevation_deg)
{
  const double ele_thresh_deg = 30; // deg

  if (elevation_deg < ele_thresh_deg) return false;
  // 視線方向が重要であることを考えて仰角最大の衛星を参照にする．
  if (elevation_deg > max_elevation_deg_)
  {
    max_elevation_deg_ = elevation_deg;
    ref_gnss_ch = count;
  }

  return true;
}

// template <typename T, size_t N>
const bool PCOEstimation::DpcoInitialEstimation(const Eigen::MatrixXd& W)
{
  static int epoch_count = 0;
  const double ddcp_res_thresh = 1e-4;

  const int N = W.rows();
  const int current_size = W_.rows();
  const int new_size = current_size + N;
  Eigen::MatrixXd W_cpy = W_;
  W_ = Eigen::MatrixXd::Zero(new_size, new_size);

  // append new matrix and observations
  // H_.block(current_size, 0, N, 3) = H;
  // V_.block(current_size, 0, N, 1) = V_Res;
  W_.topLeftCorner(current_size, current_size) = W_cpy;
  W_.block(current_size, current_size, N, N) = W;

  // std::cout << "H" << H_ << std::endl;
  // std::cout << "V" << V_ << std::endl;
  // std::cout << "W" << W_ << std::endl;

  epoch_count++;
  if (new_size >= PCO_WLS_DATA_NUM)
  {
    // Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(W_);
    // auto rank = lu_decomp.rank(); // なんかうまく計算できてなさそう．<- Wがでかすぎると行列式が0になる可能性はある．Diagonalが<1なので．
    // if (rank < W_.rows()) abort(); // 正則判定．

    Eigen::Vector3d dpco = - (H_.transpose() * W_ * H_).inverse() * (H_.transpose() * W_ * V_); // 3次元

// debug ++++++++++++++++++++++++++++++++++++++++++++++++
    Eigen::VectorXd dd_pcc_after = - H_ * dpco;
    // std::cout << "H dpco" << H_ << std::endl;
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
      dpco_mm_ = libra::Vector<3>(0);
      for (int i = 0; i < 3; i++) dpco_mm_[i] = dpco(i) * 1000;
      // 大体0.1mm以下の精度になったら収束判定をする．<- これはなぜこの値にしたのかをもう少し定量的に説明したい．PCOの推定精度は1mm以下を目指していて残差で評価するとこれくらいに相当するから，的な．
      if (pre_acc < ddcp_res_thresh && post_acc < ddcp_res_thresh)
      {
        pco_fixed_ = true;
      }
      return true;
    }

    // リセット
    V_.resize(0);
    H_.resize(0, 0);
    W_.resize(0, 0);
    epoch_count = 0;
  }
  return false;
}
