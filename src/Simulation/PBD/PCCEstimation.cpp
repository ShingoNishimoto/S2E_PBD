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

void PCCEstimation::ResizeHV(const int count)
{
  if (!pco_estimate_.GetPcoFixed()) pco_estimate_.ResizeHV(count);
  else if (!pcv_estimate_.GetPcvFixed()) pcv_estimate_.ResizeHV(count);
  else std::cout << "ERROR: something wrong" << std::endl;
}

void PCCEstimation::GetObservableInfo(const int local_pos, const int i, const int ref_j, const PBD_GnssObservation& gnss_observation, const double res_ddcp)
{
  if (!pco_estimate_.GetPcoFixed()) pco_estimate_.SetHVRaw(local_pos, i, ref_j, gnss_observation, res_ddcp);
  else if (!pcv_estimate_.GetPcvFixed()) pcv_estimate_.GetObservableInfo(local_pos, i, ref_j, gnss_observation, res_ddcp, pcc_);
  else std::cout << "ERROR: something wrong" << std::endl;
}

const bool PCCEstimation::CheckDataForEstimation(const int count, int& ref_gnss_ch, const double elevation_deg, const double r_sdcp)
{
  if (!pco_estimate_.GetPcoFixed())
  {
    if (pco_estimate_.CheckDataForEstimation(count, ref_gnss_ch, elevation_deg))
    {
      data_available_ = true;
      return true;
    }
  }
  else if (!pcv_estimate_.GetPcvFixed())
  {
    pcv_estimate_.UpdateReferenceSat(count, ref_gnss_ch, r_sdcp, elevation_deg);
    data_available_ = pcv_estimate_.data_available_;
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

// debug出力で反くて普通にログに残すとかできるように修正したい．
const bool PCCEstimation::Update(const Eigen::MatrixXd& W, const double elapsed_time)
{
  if (!pco_estimate_.GetPcoFixed())
  {
    // ここで更新する．
    if (pco_estimate_.DpcoInitialEstimation(W, elapsed_time))
    {
      pcc_->UpdatePCO(pco_estimate_.dpco_mm_);
      return true;
    }
  }
  else if (!pcv_estimate_.GetPcvFixed())
  {
    if (pcv_estimate_.Update(W, pcc_, elapsed_time))
    {
      pcc_->UpdatePCV(pcv_estimate_.dpcv_vec_mm_);
       // ステップごとに保存できるようにしたい．
      pcc_->PcvLogOutput(pcc_->out_fname_base_ + "_pcv.csv");
      pcc_->PccLogOutput(pcc_->out_fname_base_ + "_pcc.csv");
      return true;
    }
  }
  else std::cout << "ERROR: something wrong" << std::endl;

  return false;
}

