#include "PhaseCenterCorrection.hpp"
#include <Constant.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

// #define DISTANCE_BASED_INTERPOLATION
#define PIECEWISE_FUNCTION

#define PCO_WLS_DATA_NUM (1 * 11)

// PCOは左手系(north, east, up)で定義されている．左手系なのでazimuthは時計まわり．これに合うように要修正！
PhaseCenterCorrection::PhaseCenterCorrection(libra::Vector<3> pco, std::vector<double> pcv, const double azi_increment, const double ele_increment): pco_mm_(pco), pcv_mm_(pcv),
azi_increment_(azi_increment), ele_increment_(ele_increment)
{
  double azimuth;
  const int num_azi = (int)(360 / azi_increment);
  for (int i = 0; i <= num_azi; i++)
  {
    azimuth = azi_increment*i;
    azimuth_index_[azimuth] = i;
  }

  // antexファイルにはzenith angle順で並んでいるので注意．<- 天頂角に合わせた方がいいのかもしれない．
  double elevation;
  const int num_ele = (int)(90 / ele_increment);
  for (int i = 0; i <= num_ele; i++)
  {
    elevation = ele_increment*(num_ele - i);
    elevation_index_[elevation] = i;
  }

  PccLogOutput();
}

// pcoだけを指定する場合．推定用．
PhaseCenterCorrection::PhaseCenterCorrection(libra::Vector<3> pco, const double azi_increment, const double ele_increment): pco_mm_(pco),
azi_increment_(azi_increment), ele_increment_(ele_increment)
{
  const int num_azi = (int)(360 / azi_increment) + 1;
  const int num_ele = (int)(90 / ele_increment) + 1;
  pcv_mm_.assign(num_azi*num_ele, 0.0); // 0埋め
}

PhaseCenterCorrection::~PhaseCenterCorrection(){}

void PhaseCenterCorrection::PccLogOutput(void)
{
  std::ofstream ofs_("phase_center_correction.csv"); // ここはアンテナごとに指定したい．
  const int precision = 5;

  for (int i = 0; i < 3; i++)
  {
    ofs_ << std::fixed << std::setprecision(precision) << pco_mm_[i] << ",";
  }
  ofs_ << std::endl; // 改行

  const int num_ele = (int)(90 / ele_increment_) + 1;

  for (int azimuth = 0; azimuth <= 360; azimuth+=azi_increment_)
  {
    for (int elevation = 0; elevation < 90; elevation+=ele_increment_)
    {
      const int index = azimuth_index_[azimuth] * num_ele + elevation_index_[elevation];
      ofs_ << std::fixed << std::setprecision(precision) << pcv_mm_.at(index) << ",";
    }
    // 最後にカンマが入らないように調整
    const int index = azimuth_index_[azimuth] * num_ele + elevation_index_[90];
    ofs_ << std::fixed << std::setprecision(precision) << pcv_mm_.at(index) << std::endl; // 改行
  }
}

// ベクトルではなくrangeに加わる誤差としてPCCは定義される．
const double PhaseCenterCorrection::GetPCC_m(const double azimuth_deg, const double elevation_deg)
{
  // incrementの分解能で近い角度を計算
  int azi_floor =  std::floor(azimuth_deg / azi_increment_) * azi_increment_;
  int azi_ceil =  std::ceil(azimuth_deg / azi_increment_) * azi_increment_;
  int ele_floor =  std::floor(elevation_deg / ele_increment_) * ele_increment_;
  int ele_ceil =  std::ceil(elevation_deg / ele_increment_) * ele_increment_;

  // 対称点周り4点から補間して求める
  const int num_ele = 90 / ele_increment_;
  double pcv_mm_1 = pcv_mm_.at(num_ele*azimuth_index_[azi_floor] + elevation_index_[ele_floor]);
  double pcv_mm_2 = pcv_mm_.at(num_ele*azimuth_index_[azi_ceil] + elevation_index_[ele_floor]);
  double pcv_mm_3 = pcv_mm_.at(num_ele*azimuth_index_[azi_floor] + elevation_index_[ele_ceil]);
  double pcv_mm_4 = pcv_mm_.at(num_ele*azimuth_index_[azi_ceil] + elevation_index_[ele_ceil]);

#ifdef DISTANCE_BASED_INTERPOLATION
  double w_1 = 1 / sqrt(pow(azimuth_deg - azi_floor, 2.0) + pow(elevation_deg - ele_floor, 2.0));
  double w_2 = 1 / sqrt(pow(azimuth_deg - azi_ceil, 2.0) + pow(elevation_deg - ele_floor, 2.0));
  double w_3 = 1 / sqrt(pow(azimuth_deg - azi_floor, 2.0) + pow(elevation_deg - ele_ceil, 2.0));
  double w_4 = 1 / sqrt(pow(azimuth_deg - azi_ceil, 2.0) + pow(elevation_deg - ele_ceil, 2.0));
  double w_sum = (w_1 + w_2 + w_3 + w_4);
  w_1 /= w_sum; w_2 /= w_sum; w_3 /= w_sum; w_4 /= w_sum; // 正規化
#endif // DISTANCE_BASED_INTERPOLATION
#ifdef PIECEWISE_FUNCTION
  const double gamma = (azimuth_deg - azi_floor) / (azi_ceil - azi_floor);
  const double beta = (elevation_deg - ele_floor) / (ele_ceil - ele_floor);
  double w_1 = (1 - gamma) * (1 - beta);
  double w_2 = gamma * (1 - beta);
  double w_3 = (1 - gamma) * beta;
  double w_4 = gamma * beta;
#endif // PIECEWISE_FUNCTION

  // 距離に応じた重み付き平均 <- 角度に対してユークリッド距離を定義するのは微妙な気がするので，論文の手法に従う方がよさそう．
  double target_pcv = w_1*pcv_mm_1 + w_2*pcv_mm_2 + w_3*pcv_mm_3 + w_4*pcv_mm_4;

  const double azi_rad = azimuth_deg * libra::deg_to_rad;
  const double ele_rad = elevation_deg * libra::deg_to_rad;
  std::vector<double> e_vec = { cos(ele_rad) * cos(azi_rad),
                                cos(ele_rad) * sin(azi_rad),
                                sin(ele_rad) };
  // この時PCOはARP固定座標系であるが，Azimuth，Elevationもコンポ固定の座標系であるため特に変換を入れてない．コンポ座標系とARP固定座標が一致してるかどうかは要注意．
  double pcc;
  pcc = -(pco_mm_[0]*e_vec.at(0) + pco_mm_[1]*e_vec.at(1) + pco_mm_[2]*e_vec.at(2)) + target_pcv;
  pcc /= 1000; // mに変換
  return pcc;
}

// ここはPCOEstimationクラスとして別だしした方がいいかも
// template <typename T, size_t N>
void PhaseCenterCorrection::DpcoInitialEstimation(const Eigen::MatrixXd& H, const Eigen::VectorXd& V_Res, const Eigen::MatrixXd& W)
{
  static int epoch_count = 0;
  const double ddcp_res_thresh = 1e-4;

  const int N = V_Res.rows();
  const int current_size = V_Res_dpco_.rows();
  const int new_size = current_size + N;
  H_dpco_.conservativeResize(new_size, 3);
  V_Res_dpco_.conservativeResize(new_size);
  Eigen::MatrixXd W_cpy = W_dpco_;
  W_dpco_ = Eigen::MatrixXd::Zero(new_size, new_size);

  // append new matrix and observations
  H_dpco_.block(current_size, 0, N, 3) = H;
  V_Res_dpco_.block(current_size, 0, N, 1) = V_Res;
  W_dpco_.topLeftCorner(current_size, current_size) = W_cpy;
  W_dpco_.block(current_size, current_size, N, N) = W;

  // std::cout << "H" << H_dpco_ << std::endl;
  // std::cout << "V" << V_Res_dpco_ << std::endl;
  // std::cout << "W" << W_dpco_ << std::endl;

  epoch_count++;
  if (new_size >= PCO_WLS_DATA_NUM)
  {
    Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(W_dpco_);
    auto rank = lu_decomp.rank(); // なんかうまく計算できてなさそう．<- Wがでかすぎると行列式が0になる可能性はある．Diagonalが<1なので．
    // if (rank < W_dpco_.rows()) abort(); // 正則判定．

    Eigen::Vector3d dpco = - (H_dpco_.transpose() * W_dpco_ * H_dpco_).inverse() * (H_dpco_.transpose() * W_dpco_ * V_Res_dpco_); // 3次元

// debug ++++++++++++++++++++++++++++++++++++++++++++++++
    Eigen::VectorXd dd_pcc_after = - H_dpco_ * dpco;
    // std::cout << "H dpco" << H_dpco_ << std::endl;
    // std::cout << "after estimated ddpcc" << dd_pcc_after << std::endl;
    // std::cout << "raw ddpcc" << V_Res_dpco_ << std::endl;
    double pre_acc = 0;
    for (int i = 0; i < V_Res_dpco_.rows(); i++) pre_acc += pow(V_Res_dpco_(i), 2.0);
    pre_acc = sqrt(pre_acc) / V_Res_dpco_.rows(); // RMS
    // std::cout << "pre accuracy" << pre_acc << std::endl;

    Eigen::VectorXd post_res = V_Res_dpco_ - dd_pcc_after;
    double post_acc = 0;
    for (int i = 0; i < post_res.rows(); i++) post_acc += pow(post_res(i), 2.0);
    post_acc = sqrt(post_acc) / post_res.rows();
    // std::cout << "post accuracy" << post_acc << std::endl;
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // 改善した時だけ更新する．
    if (pre_acc > post_acc)
    {
      libra::Vector<3> dpco_mm(0);
      for (int i = 0; i < 3; i++) dpco_mm[i] = dpco(i) * 1000;
      UpdatePCO(dpco_mm);
      // 大体0.1mm以下の精度になったら収束判定をする．
      if (pre_acc < ddcp_res_thresh && post_acc < ddcp_res_thresh)
      {
        pco_fixed_ = true;
      }
    }

    // リセット
    V_Res_dpco_.resize(0);
    H_dpco_.resize(0, 0);
    W_dpco_.resize(0, 0);
    epoch_count = 0;
  }
}
