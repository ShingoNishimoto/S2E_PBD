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

  void ResizeHV(const int count);
  void GetObservableInfo(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp);
  const bool CheckDataForEstimation(const int count, int& ref_gnss_ch, const double elevation_deg, const double r_sdcp);
  inline const bool GetEstimationFinish(void) const { return estimation_finish_; }
  void InitializeRefInfo(void);
  const bool Update(const Eigen::MatrixXd& W); // const Eigen::VectorXd& V,

  bool data_available_ = false;

 protected:
  PhaseCenterCorrection* pcc_;
  PCOEstimation pco_estimate_;
  PCVEstimation pcv_estimate_;

  bool estimation_finish_ = false;
};