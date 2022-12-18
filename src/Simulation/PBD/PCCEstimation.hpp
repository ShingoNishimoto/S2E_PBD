#pragma once

#include <map>
#include "PhaseCenterCorrection.hpp"
#include "PCOEstimation.hpp"
#include "PCVEstimation.hpp"

using libra::Matrix;
using libra::Vector;
using std::vector;

class PCCEstimation
{
 public:
  PCCEstimation(PhaseCenterCorrection* pcc, const std::string fname);
  PCCEstimation();
  ~PCCEstimation();

  void ResizeH(const int count);
  void SetHRaw(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation);
  const bool CheckDataForEstimation(const int count, int& ref_gnss_ch, const double elevation_deg, const double r_sdcp);
  inline const bool GetEstimationFinish(void) const { return estimation_finish_; }
  void InitializeRefInfo(void);
  void Update(const Eigen::VectorXd& V, const Eigen::MatrixXd& W);

 protected:
  PhaseCenterCorrection* pcc_;
  PCOEstimation pco_estimate_;
  PCVEstimation pcv_estimate_;

  bool estimation_finish_ = false;
};
