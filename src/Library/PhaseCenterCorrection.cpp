#include "PhaseCenterCorrection.hpp"
#include <Constant.hpp>

PhaseCenterCorrection::PhaseCenterCorrection(libra::Vector<3> pco, std::vector<double> pcv, const double azi_increment, const double ele_increment): pco_(pco), pcv_(pcv),
azi_increment_(azi_increment), ele_increment_(ele_increment)
{
  double azimuth;
  const int num_azi = (int)(360 / azi_increment);
  for (int i = 0; i <= num_azi; i++)
  {
    azimuth = azi_increment*i;
    azimuth_index_[azimuth] = i;
  }

  // antexファイルにはzenith angle順で並んでいるので注意．
  double elevation;
  const int num_ele = (int)(90 / ele_increment);
  for (int i = 0; i <= num_ele; i++)
  {
    elevation = ele_increment*(num_ele - i);
    elevation_index_[elevation] = i;
  }
}

const libra::Vector<3> PhaseCenterCorrection::GetPCC(const double azimuth_deg, const double elevation_deg)
{
  // incrementの分解能で近い角度を計算
  int azi_floor =  std::floor(azimuth_deg / azi_increment_) * azi_increment_;
  int azi_ceil =  std::ceil(azimuth_deg / azi_increment_) * azi_increment_;
  int ele_floor =  std::floor(elevation_deg / ele_increment_) * ele_increment_;
  int ele_ceil =  std::ceil(elevation_deg / ele_increment_) * ele_increment_;

  // 対称点周り4点から保管して求める
  const int num_ele = 90 / ele_increment_;
  double pcv_1 = pcv_.at(num_ele*azimuth_index_[azi_floor] + elevation_index_[ele_floor]);
  double pcv_2 = pcv_.at(num_ele*azimuth_index_[azi_ceil] + elevation_index_[ele_floor]);
  double pcv_3 = pcv_.at(num_ele*azimuth_index_[azi_floor] + elevation_index_[ele_ceil]);
  double pcv_4 = pcv_.at(num_ele*azimuth_index_[azi_ceil] + elevation_index_[ele_ceil]);

  double w_1 = 1 / sqrt(pow(azimuth_deg - azi_floor, 2.0) + pow(elevation_deg - ele_floor, 2.0));
  double w_2 = 1 / sqrt(pow(azimuth_deg - azi_ceil, 2.0) + pow(elevation_deg - ele_floor, 2.0));
  double w_3 = 1 / sqrt(pow(azimuth_deg - azi_floor, 2.0) + pow(elevation_deg - ele_ceil, 2.0));
  double w_4 = 1 / sqrt(pow(azimuth_deg - azi_ceil, 2.0) + pow(elevation_deg - ele_ceil, 2.0));

  // 距離に応じた重み付き平均
  double target_pcv = (w_1*pcv_1 + w_2*pcv_2 + w_3*pcv_3 + w_4*pcv_4) / (w_1 + w_2 + w_3 + w_4);

  libra::Vector<3> pcc(0);
  pcc[0] = pco_[0] + target_pcv * cos(elevation_deg * libra::deg_to_rad) * cos(azimuth_deg * libra::deg_to_rad);
  pcc[1] = pco_[1] + target_pcv * cos(elevation_deg * libra::deg_to_rad) * sin(azimuth_deg * libra::deg_to_rad);
  pcc[2] = pco_[2] + target_pcv * sin(elevation_deg * libra::deg_to_rad);
  for (int i = 0; i < 3; i++) pcc[i] /= 1000; // mに変換
  return pcc;
}
