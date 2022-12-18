#pragma once

#include "../../Library/VectorTool.hpp"
#include <Matrix.hpp>
#include <map>
#include "PBD_GnssObservation.h"

enum class PCV_METHOD
{
  SPHERE,   // Direct approach, use spherical harmonics
  ZERNIKE,  // Direct approach, use Zernike polynomial
  RESIDUAL, // Residual approach
};

using libra::Matrix;
using libra::Vector;
using std::vector;

class PCVEstimation
{
 public:
  PCVEstimation(const std::string fname);
  PCVEstimation();
  ~PCVEstimation();

  void SetHRaw(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation);
  void UpdateReferenceSat(const int count, int& ref_gnss_ch, const double r_sdcp);
  inline const bool GetPcvFixed(void) const { return pcv_fixed_; }
  const bool Update(const Eigen::VectorXd& V, const Eigen::MatrixXd& W, const double azi_increment, const double ele_increment);
  inline void ResizeH(const int count) { H_.conservativeResize((int)H_.rows() + count - 1 , (degree_ + 2) * (degree_ + 1)); }

  double min_variance_ = 1e18; // 初期値は適当に大きな値．
  vector<double> dpcv_vec_mm_;

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
  Eigen::VectorXd CS_vec_;
  Eigen::VectorXd V_;
  Eigen::MatrixXd H_;
  Eigen::MatrixXd W_;

  void SphericalHarmonicsInitialization(const std::string fname);
  const Eigen::MatrixXd CalcLegendreCoeff(const double azi_rad, const double ele_rad);
  const bool WeightedLeastSquare(const Eigen::VectorXd& V, const Eigen::MatrixXd& W, const double azi_increment, const double ele_increment);
  void SetPcvVecFromSHModel(const double azi_increment, const double ele_increment);
  void RemoveZeroCols(Eigen::MatrixXd& H);
};
