#include "PCCEstimation.hpp"
#include <Constant.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Interface/InitInput/IniAccess.h"

PCCEstimation::PCCEstimation(PhaseCenterCorrection* pcc, const std::string fname): pcc_(pcc)
{
  pco_estimate_ = PCOEstimation();
  pcv_estimate_ = PCVEstimation(fname);
}

PCCEstimation::PCCEstimation(){}

PCCEstimation::~PCCEstimation(){}

void PCCEstimation::ResizeH(const int count)
{
  if (!pco_estimate_.GetPcoFixed()) pco_estimate_.ResizeH(count);
  else if (!pcv_estimate_.GetPcvFixed()) pcv_estimate_.ResizeH(count);
  else std::cout << "ERROR: something wrong" << std::endl;
}

void PCCEstimation::SetHRaw(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation)
{
  if (!pco_estimate_.GetPcoFixed()) pco_estimate_.SetHRaw(local_pos, i, ref_j, gnss_observation);
  else if (!pcv_estimate_.GetPcvFixed()) pcv_estimate_.SetHRaw(local_pos, i, ref_j, gnss_observation);
  else std::cout << "ERROR: something wrong" << std::endl;
}

const bool PCCEstimation::CheckDataForEstimation(const int count, int& ref_gnss_ch, const double elevation_deg, const double r_sdcp)
{
  if (!pco_estimate_.GetPcoFixed())
  {
    if (pco_estimate_.CheckDataForEstimation(count, ref_gnss_ch, elevation_deg)) return true;
  }
  else if (!pcv_estimate_.GetPcvFixed())
  {
    pcv_estimate_.UpdateReferenceSat(count, ref_gnss_ch, r_sdcp, elevation_deg);
    return true;
  }
  else std::cout << "ERROR: something wrong" << std::endl;

  return false;
}

void PCCEstimation::InitializeRefInfo(void)
{
  pco_estimate_.max_elevation_deg_ = 0.0;
  pcv_estimate_.max_elevation_deg_ = 0.0;
  pcv_estimate_.min_variance_ = 1e18;
}

const bool PCCEstimation::Update(const Eigen::VectorXd& V, const Eigen::MatrixXd& W)
{
  if (!pco_estimate_.GetPcoFixed())
  {
    // ここで更新する．
    if (pco_estimate_.DpcoInitialEstimation(V, W))
    {
      pcc_->UpdatePCO(pco_estimate_.dpco_mm_);
      return true;
    }
  }
  else if (!pcv_estimate_.GetPcvFixed())
  {
    if (pcv_estimate_.Update(V, W, pcc_->azi_increment_, pcc_->ele_increment_))
    {
      pcc_->UpdatePCV(pcv_estimate_.dpcv_vec_mm_);
      pcc_->PccLogOutput(pcc_->out_fname_base_ + ".csv"); // ステップごとに保存できるようにしたい．
      return true;
    }
  }
  else std::cout << "ERROR: something wrong" << std::endl;

  return false;
}

