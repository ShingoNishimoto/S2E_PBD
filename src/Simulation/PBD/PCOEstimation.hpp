#pragma once

#include "../../Library/VectorTool.hpp"
#include <Matrix.hpp>
#include <map>
#include "PBD_GnssObservation.h"

using libra::Matrix;
using libra::Vector;
using std::vector;

class PCOEstimation
{
 public:
  PCOEstimation();
  ~PCOEstimation();

  // template <typename T, size_t N>
  const bool DpcoInitialEstimation(const Eigen::VectorXd& V_Res, const Eigen::MatrixXd& W);
  const bool GetPcoFixed(void) const { return pco_fixed_; }
  void SetHRaw(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation);
  inline void ResizeH(const int count) { H_.conservativeResize((int)H_.rows() + count - 1 , 3); }
  const bool CheckDataForEstimation(const int count, int& ref_gnss_ch, const double elevation_deg);

  double max_elevation_deg_ = 0.0;
  libra::Vector<3> dpco_mm_;

 protected:

  Eigen::VectorXd V_;
  Eigen::MatrixXd H_;
  Eigen::MatrixXd W_;
  bool pco_fixed_ = false;

  const Eigen::VectorXd WeightedLeastSquare(const Eigen::MatrixXd& H, const Eigen::VectorXd& V, const Eigen::MatrixXd& W);
};
