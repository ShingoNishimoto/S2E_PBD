#include "PCVEstimation.hpp"
#include <Constant.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Interface/InitInput/IniAccess.h"

PCVEstimation::PCVEstimation(PCV_METHOD method, const std::string fname): method_(method)
{
  switch (method)
  {
  case PCV_METHOD::SPHERE:
    SphericalHarmonicsInitialization(fname);
    break;

  case PCV_METHOD::ZERNIKE:
    // TODO
    break;

  case PCV_METHOD::RESIDUAL:
    // TODO
    break;

  default:
    break;
  }
}

void PCVEstimation::SphericalHarmonicsInitialization(const std::string fname)
{
  IniAccess pcv_conf(fname);
  degree_ = pcv_conf.ReadInt("SPHERE", "degree");
  // coefficients
  c_.assign(degree_ + 1, vector<double>(degree_ + 1, 0.0));
  s_.assign(degree_ + 1, vector<double>(degree_ + 1, 0.0));
  CS_vec_ = Eigen::VectorXd::Zero((degree_ + 1) * (degree_ + 1));
}

void PCVEstimation::CalcLegendreCoeff(const double azi_rad, const double ele_rad)
{
  x_ = cos(ele_rad) * cos(azi_rad);
  y_ = cos(ele_rad) * sin(azi_rad);
  z_ = sin(ele_rad);

  // Calc V and W
  int degree_vw = degree_ + 1;
  vector<vector<double>> v(degree_vw + 1, vector<double>(degree_vw + 1, 0.0));
  vector<vector<double>> w(degree_vw + 1, vector<double>(degree_vw + 1, 0.0));
  // n_=m_=0
  v[0][0] = 1.0;
  w[0][0] = 0.0;
  m_ = 0;

  Eigen::MatrixXd VW = Eigen::MatrixXd::Zero(1, (degree_ + 1) * (degree_ + 1)); // 推定に使う方はCSと同じ次元．

  while (m_ < degree_vw) {
    for (n_ = m_ + 1; n_ <= degree_vw; n_++) {
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
  // このV, Wをもとに係数ベクトルを作成し，それをもとに推定していくことになる．計算量やばそう．．．低次数順に並べていけばシンプルにできそう．
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
