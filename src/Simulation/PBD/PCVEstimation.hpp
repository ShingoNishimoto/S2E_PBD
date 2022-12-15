#pragma once

#include "../../Library/VectorTool.hpp"
#include <Matrix.hpp>
#include <map>

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
  PCVEstimation(PCV_METHOD method, const std::string fname);
  ~PCVEstimation();

 protected:
  PCV_METHOD method_;

  // SH
  int degree_;
  vector<vector<double>> c_;
  vector<vector<double>> s_;
  int n_ = 0, m_ = 0;
  double x_ = 0.0, y_ = 0.0, z_ = 0.0;
  void v_w_nn_update(double *v_nn, double *w_nn, const double v_prev, const double w_prev);
  void v_w_nm_update(double *v_nm, double *w_nm, const double v_prev, const double w_prev, const double v_prev2, const double w_prev2);

  Eigen::VectorXd CS_vec_;
  Eigen::VectorXd V_;
  Eigen::MatrixXd H_;
  Eigen::MatrixXd W_;

  void SphericalHarmonicsInitialization(const std::string fname);
  void CalcLegendreCoeff(const double azi_rad, const double ele_rad);

};
