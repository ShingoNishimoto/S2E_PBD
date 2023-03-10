#include "PCVEstimation.hpp"
#include <Constant.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include "Interface/InitInput/IniAccess.h"

#define SH_WLS_DATA_MULTIPLE (40)
#define RES_MEAN_DATA_NUM (5)

PCV_GnssDirection::PCV_GnssDirection(const int azimuth, const int elevation): azimuth_(azimuth), elevation_(elevation)
{}

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
    ResidualInitialization(fname);
  }
  else
  {
    std::cout << "ERROR: PCV estimation method error!" << std::endl;
  }
}

PCVEstimation::PCVEstimation(){}

PCVEstimation::~PCVEstimation(){}


void PCVEstimation::GetObservableInfo(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp, PhaseCenterCorrection* pcc)
{
  switch (method_)
  {
  case PCV_METHOD::SPHERE:
    SetHVRaw(local_pos, i, ref_j, gnss_observation, res_ddcp);
    break;

  case PCV_METHOD::ZERNIKE:
    // TODO
    break;

  case PCV_METHOD::RESIDUAL:
    SetGnssInfo(local_pos, i, ref_j, gnss_observation, res_ddcp, pcc);
    break;

  default:
    // ERROR
    break;
  }
}

// const bool PCVEstimation::CheckDataForEstimation(const int count, int& ref_gnss_ch, const double r_sdcp, const double elevation_deg)
// {
//   switch (method_)
//   {
//   case PCV_METHOD::SPHERE:
//     UpdateReferenceSat(count, ref_gnss_ch, r_sdcp, elevation_deg);
//     return true;

//   case PCV_METHOD::ZERNIKE:
//     // TODO
//     return false;

//   case PCV_METHOD::RESIDUAL:
//     UpdateReferenceSat(count, ref_gnss_ch, r_sdcp, elevation_deg);
//     return true;

//   default:
//     return false;
//   }
// }

void PCVEstimation::UpdateReferenceSat(const int count, int& ref_gnss_ch, const double r_sdcp, const double elevation_deg)
{
  // // 観測誤差が小さな衛星を参照衛星とする．
  // if (r_sdcp < min_variance_)
  // {
  //   min_variance_ = r_sdcp;
  //   ref_gnss_ch = count;
  // }

  // 仰角を参照するパターン
  if (elevation_deg > max_elevation_deg_)
  {
    max_elevation_deg_ = elevation_deg;
    ref_gnss_ch = count;
  }

  data_available_ = true;
  if (method_ == PCV_METHOD::RESIDUAL)
  {
    // 最大仰角がTBD以下の時は使えない．
    if (max_elevation_deg_ < 75) data_available_ = false;
  }
}

// 更新が入った時だけtrueを返す．
const bool PCVEstimation::Update(const Eigen::MatrixXd& W, PhaseCenterCorrection* pcc, const double elapsed_time)
{
  bool result = false;
  switch (method_)
  {
  case PCV_METHOD::SPHERE:
    // return WeightedLeastSquare(W, azi_increment, ele_increment);
    // Wを実際のものにするとうまく行かないので一旦単位行列で実施する．
    result = WeightedLeastSquare(Eigen::MatrixXd::Identity(W.rows(), W.cols()), pcc->azi_increment_, pcc->ele_increment_);
    break;

  case PCV_METHOD::ZERNIKE:
    // TODO
    break;

  case PCV_METHOD::RESIDUAL:
    result = ResidualBasedUpdate(W, pcc);
    break;

  default:
    return false;
  }

  if (result)
  {
    std::cout << "PCV updated! at " << elapsed_time << std::endl;
  }
  return result;
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

void PCVEstimation::SetHVRaw(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp)
{
  const double azi_rad_i = gnss_observation.GetGnssAzimuthDeg(i) * libra::deg_to_rad;
  const double ele_rad_i = gnss_observation.GetGnssElevationDeg(i) * libra::deg_to_rad;
  const double azi_rad_j = gnss_observation.GetGnssAzimuthDeg(ref_j) * libra::deg_to_rad;
  const double ele_rad_j = gnss_observation.GetGnssElevationDeg(ref_j) * libra::deg_to_rad;
  const Eigen::MatrixXd Hi = CalcLegendreCoeff(azi_rad_i, ele_rad_i) - CalcLegendreCoeff(azi_rad_j, ele_rad_j);
  const int pos = local_pos + W_.rows();
  H_.block(pos, 0, 1, Hi.cols()) = Hi;
  V_(pos) = res_ddcp;
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
  V_ = Eigen::VectorXd::Zero(1);
  W_ = Eigen::MatrixXd::Identity(1, 1);// * 1e18; // constraintなので不確定性は0として．

  if (method_ == PCV_METHOD::SPHERE)
  {
    const int coef_size = (degree_ + 2) * (degree_ + 1);
    H_ = Eigen::MatrixXd::Zero(1, coef_size);

    vector<double> pn0(degree_ + 1);
    pn0[0] = 1.0; pn0[1] = 1.0;
    for (n_ = 0; n_ <= degree_; n_++)
    {
      if (n_ >= 2) p_n_0_update(&pn0[n_], pn0[n_ - 1], pn0[n_ - 2]);
      H_(0, n_ * (n_ + 1)) = pn0[n_];
    }
    // std::cout << "H" << H_ << std::endl;
  }
}

// TODO: マスク角によって見えてない部分は観測量が得られないので，その影響で誤差が大きくなっている可能性がある．ここはモデルに入れないようにすることで精度を上げることはできる可能性があるので修正する．
// jとiの差分なので，衛星のLOSベクトルの差の角度が大きい組み合わせじゃないと実施しないとかをやらないと精度良くするのはムズイのかも？実は組み合わせやからjを一つに固定する必要はないな．
const bool PCVEstimation::WeightedLeastSquare(const Eigen::MatrixXd& W, const double azi_increment, const double ele_increment)
{
  const int N = W.rows();
  const int current_size = W_.rows();
  const int new_size = current_size + N;
  const int coef_num = (degree_ + 2) * (degree_ + 1);
  H_.conservativeResize(new_size, coef_num);
  Eigen::MatrixXd W_cpy = W_;
  W_ = Eigen::MatrixXd::Zero(new_size, new_size);

  // append new matrix and observations
  // H_.block(current_size, 0, N, 3) = H;
  // V_.block(current_size, 0, N, 1) = V;
  W_.topLeftCorner(current_size, current_size) = W_cpy;
  W_.block(current_size, current_size, N, N) = W;

  // std::cout << "H" << H_ << std::endl;
  // std::cout << "V" << V_ << std::endl;
  // std::cout << "W" << W_ << std::endl;

  static double ddcp_res_thresh = 1e-4;
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
      // 収束判定の閾値はどうする？
      if (post_acc < ddcp_res_thresh) // pre_acc < ddcp_res_thresh &&
      {
        pcv_fixed_ = true;
        ddcp_res_thresh *= 0.8;
      }
      InitializeVHW();
      return true;
    }

    // リセット
    InitializeVHW();
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
      double elevation = ele_increment * (num_ele - j);
      Eigen::MatrixXd VW = CalcLegendreCoeff(azimuth * libra::deg_to_rad, elevation * libra::deg_to_rad);
      RemoveZeroCols(VW);
      const double dpcv_mm = (VW * CS_vec_)(0) * 1000.0; // mmに変換．
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

void PCVEstimation::ResidualInitialization(const std::string fname)
{
  IniAccess pcv_conf(fname);
  res_azi_increment_ = pcv_conf.ReadDouble("RESIDUAL", "azi_increment");
  res_ele_increment_ = pcv_conf.ReadDouble("RESIDUAL", "ele_increment");

  const int num_azi = (int)(360 / res_azi_increment_);
  const int num_ele = (int)(90 / res_ele_increment_);
  // groupをelevation方向だけにした方がいいかも．そうすればもう少し増えそう．
  res_mm_vec_ = vector<vector<double>>((num_azi + 1) * (num_ele + 1));
  weight_vec_ = vector<vector<double>>((num_azi + 1) * (num_ele + 1));

  for (int i = 0; i <= num_azi; i++)
  {
    double azimuth = res_azi_increment_*i;
    res_azimuth_index_[azimuth] = i;
  }

  for (int i = 0; i <= num_ele; i++)
  {
    double elevation = res_ele_increment_*(num_ele - i);
    res_elevation_index_[elevation] = i;
  }

  InitializeVHW();
}

// FIXME: ここは来たまんまで埋めていくと偏りが出るので工夫が必要．時間とのトレードオフではあるが．
void PCVEstimation::SetGnssInfo(const int ch, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp, PhaseCenterCorrection* pcc)
{
  const double ref_azimuth = gnss_observation.GetGnssAzimuthDeg(ref_j);
  const double ref_elevation = gnss_observation.GetGnssElevationDeg(ref_j);
  const double ref_pcv_mm = pcc->GetPCV_mm(ref_azimuth, ref_elevation);

  // roundを使って近いグリッド点の値とみなす．
  const double azimuth = gnss_observation.GetGnssAzimuthDeg(i);
  const double elevation = gnss_observation.GetGnssElevationDeg(i);
  int round_azi = std::round(azimuth / res_azi_increment_) * res_azi_increment_; if (round_azi >= 360) round_azi = 0;
  const int round_ele = std::round(elevation/ res_ele_increment_) * res_ele_increment_;
  const int num_ele = (int)(90 / res_ele_increment_) + 1;

  // std::cout << "azi: " << azimuth << ", ele: " << elevation << ", DDCP residual: " << res_ddcp * 1000 << std::endl;

  // 参照部分のPCVとの和をPCV観測量とみなして追加．
  // 天頂方向の残差はスキップするようにすべきなきがする．
  res_mm_vec_.at(num_ele*res_azimuth_index_.at(round_azi) + res_elevation_index_.at(round_ele)).push_back(res_ddcp * 1000 + ref_pcv_mm);

  PCV_GnssDirection dir = PCV_GnssDirection(round_azi, round_ele);
  observable_info_vec_[dir].push_back(ch);
}

const bool PCVEstimation::ResidualBasedUpdate(const Eigen::MatrixXd& W, PhaseCenterCorrection* pcc)
{
  const int pcv_num_azi = (int)(360 / pcc->azi_increment_) + 1;
  const int pcv_num_ele = (int)(90 / pcc->ele_increment_) + 1;
  dpcv_vec_mm_.assign((pcv_num_azi) * (pcv_num_ele), 0.0); // 初めに0埋め

  const int num_ele = (int)(90 / res_ele_increment_) + 1;
  bool updated = false;
  for (const auto& info : observable_info_vec_)
  {
    const PCV_GnssDirection& dir = info.first;
    const int index = num_ele*res_azimuth_index_.at(dir.azimuth_) + res_elevation_index_.at(dir.elevation_);
    for (const int& ch : info.second)
    {
      weight_vec_.at(index).push_back(W(ch, ch));
    }
    if (res_mm_vec_.at(index).size() != weight_vec_.at(index).size()) abort();

    // 一定数たまれば平均をとる．
    if (res_mm_vec_.at(index).size() >= RES_MEAN_DATA_NUM)
    {
      double dpcv_mm = CalcAverageDDCPResidual(index);

      // 複数gridに該当する場合は以下のように更新する．
      int azi_floor =  std::round((dir.azimuth_ - res_azi_increment_ * 0.5) / pcc->azi_increment_) * pcc->azi_increment_;
      int azi_ceil =  std::round((dir.azimuth_ + res_azi_increment_ * 0.5 - pcc->azi_increment_) / pcc->azi_increment_) * pcc->azi_increment_;
      int ele_floor =  std::round((dir.elevation_ - res_ele_increment_ * 0.5) / pcc->ele_increment_) * pcc->ele_increment_; if (ele_floor < 0) ele_floor = 0;
      int ele_ceil =  std::round((dir.elevation_ + res_ele_increment_ * 0.5 - pcc->ele_increment_) / pcc->ele_increment_) * pcc->ele_increment_; if (ele_ceil > 90) ele_ceil = 90;

      for (int azimuth = azi_floor; azimuth <= azi_ceil; azimuth+=pcc->azi_increment_)
      {
        for (int elevation = ele_floor; elevation <= ele_ceil; elevation+=pcc->ele_increment_)
        {
          // ここで角度の値域を直す．
          int azi_fixed = azimuth;
          if (azimuth < 0) azi_fixed += 360;
          else if (azimuth > 360) azi_fixed -= 360;
          const int dpcv_index = pcc->azimuth_index_.at(azi_fixed) * pcv_num_ele + pcc->elevation_index_.at(elevation);
          dpcv_vec_mm_.at(dpcv_index) = dpcv_mm;
          // 360にも同一の値を追加．
          if (azi_fixed == 0) dpcv_vec_mm_.at(pcc->azimuth_index_.at(360) * pcv_num_ele + pcc->elevation_index_.at(elevation)) = dpcv_mm;
        }
      }
      // クリア
      res_mm_vec_.at(index).clear(); weight_vec_.at(index).clear();
      updated = true;
    }
  }

  observable_info_vec_.clear();
  // 更新
  // pcc->UpdatePCV(dpcv_vec_mm_);
  return updated;
}

const double PCVEstimation::CalcAverageDDCPResidual(const int index)
{
  bool check_finish = false;
  vector<double> residuals = res_mm_vec_.at(index);
  vector<double> weights = weight_vec_.at(index);

  double dpcv_mm;
  while(!check_finish)
  {
    dpcv_mm = 0;
    double weight_sum = 0;
    const int data_num = residuals.size();
    if (data_num < 3) return 0; // 外れ値が多い場合は0を返す．

    for (int i = 0; i < data_num; i++)
    {
      dpcv_mm += residuals.at(i) * weights.at(i);
      weight_sum += weights.at(i);
    }
    dpcv_mm /= weight_sum; // mean

    double variance = 0;
    for (int i = 0; i < data_num; i++)
    {
      variance += pow(residuals.at(i) - dpcv_mm, 2.0);
    }
    variance /= data_num;
    const double sigma = sqrt(variance);
    check_finish = true;

    // 極端に大きな値をはじきたい．
    for (auto it = residuals.begin(); it != residuals.end();)
    {
      // 2 sigma区間に入っていれば正常とする．
      if (*it > dpcv_mm + 2 * sigma || *it < dpcv_mm - 2 * sigma)
      {
        it = residuals.erase(it);
        size_t i = std::distance(residuals.begin(), it);
        weights.erase(weights.begin() + i); // weightに関しても削除
        check_finish = false;
      }
      else ++it;
    }
  }

  if (!std::isfinite(dpcv_mm))
  {
    std::cout << "pcv error!" << std::endl;
    // abort();
    return 0;
  }
  return dpcv_mm;
}
