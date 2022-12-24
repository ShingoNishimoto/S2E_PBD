#pragma once

#include "../../Library/VectorTool.hpp"
#include <Matrix.hpp>
#include <map>
#include "PBD_GnssObservation.h"
#include "PhaseCenterCorrection.hpp"

enum class PCV_METHOD
{
  SPHERE,   // Direct approach, use spherical harmonics
  ZERNIKE,  // Direct approach, use Zernike polynomial
  RESIDUAL, // Residual approach
};

using libra::Matrix;
using libra::Vector;
using std::vector;

class PCV_GnssDirection
{
 public:
  PCV_GnssDirection(const int azimuth, const int elevation);
  int azimuth_;
  int elevation_;

  // mapのkeyとして指定するために定義
  bool operator<(const PCV_GnssDirection& rhs) const
  {
    return elevation_ < rhs.elevation_;
  }
};

class PCVEstimation
{
 public:
  PCVEstimation(const std::string fname);
  PCVEstimation();
  ~PCVEstimation();

  void GetObservableInfo(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp, PhaseCenterCorrection* pcc);
  // const bool CheckDataForEstimation(const int count, int& ref_gnss_ch, const double r_sdcp, const double elevation_deg);
  void UpdateReferenceSat(const int count, int& ref_gnss_ch, const double r_sdcp, const double elevation_deg);
  inline const bool GetPcvFixed(void) const { return pcv_fixed_; }
  const bool Update(const Eigen::MatrixXd& W, PhaseCenterCorrection* pcc, const double elapsed_time);
  inline void ResizeHV(const int count) { if (method_ == PCV_METHOD::SPHERE) { H_.conservativeResize((int)H_.rows() + count - 1 , (degree_ + 2) * (degree_ + 1)); V_.conservativeResize((int)V_.rows() + count - 1); } } // この辺の処理は汎化したい．

  double min_variance_ = 1e18; // 初期値は適当に大きな値．
  double max_elevation_deg_ = 0.0;
  vector<double> dpcv_vec_mm_; // PCVと同様のindex順．
  bool data_available_ = false;

 protected:
  PCV_METHOD method_;
  bool pcv_fixed_ = false;

  // SH
  int degree_;
  int wsl_data_num_; // WSLの時に使う観測データ数．
  vector<vector<double>> c_;
  vector<vector<double>> s_;
  int n_ = 0, m_ = 0;
  double x_ = 0.0, y_ = 0.0, z_ = 0.0;
  void v_w_nn_update(double *v_nn, double *w_nn, const double v_prev, const double w_prev);
  void v_w_nm_update(double *v_nm, double *w_nm, const double v_prev, const double w_prev, const double v_prev2, const double w_prev2);
  void p_n_0_update(double *p_n0, const double p_prev, const double p_prev2);
  void InitializeVHW(void);
  void SetHVRaw(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp);
  Eigen::VectorXd CS_vec_;
  Eigen::VectorXd V_;
  Eigen::MatrixXd H_;
  Eigen::MatrixXd W_;

  void SphericalHarmonicsInitialization(const std::string fname);
  const Eigen::MatrixXd CalcLegendreCoeff(const double azi_rad, const double ele_rad);
  const bool WeightedLeastSquare(const Eigen::MatrixXd& W, const double azi_increment, const double ele_increment);
  void SetPcvVecFromSHModel(const double azi_increment, const double ele_increment);
  void RemoveZeroCols(Eigen::MatrixXd& H);

  // Residual Approach

  vector<vector<double>> res_mm_vec_;
  std::map<PCV_GnssDirection, vector<int>> observable_info_vec_;
  vector<vector<double>> weight_vec_;
  double res_azi_increment_; // pccのモデルと違っても対応できるようにしている．
  double res_ele_increment_; // pccのモデルと違っても対応できるようにしている．
  std::map<int, int> res_azimuth_index_;   // azimuth角を入力としたindex取得用辞書
  std::map<int, int> res_elevation_index_; // elevation角を入力としたindex取得用辞書
  double azimuth_, elevation_;
  double ref_azimuth_, ref_elevation_;
  void ResidualInitialization(const std::string fname);
  void SetGnssInfo(const int ch, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp, PhaseCenterCorrection* pcc);
  const bool ResidualBasedUpdate(const Eigen::MatrixXd& W, PhaseCenterCorrection* pcc);
  const double CalcAverageDDCPResidual(const int index);
};
