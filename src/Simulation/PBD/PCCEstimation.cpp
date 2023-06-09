#include "PCCEstimation.hpp"
#include <Constant.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "Interface/InitInput/IniAccess.h"

// #define ITERATION
// #define WITHOUT_HORIZON_UPDATE

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
#ifdef WITHOUT_HORIZON_UPDATE
      // 水平構成分は更新しないようにする．<- ただこれだと他成分の影響が入ってる気もする．
      pco_estimate_.dpco_mm_[0] = 0.0;
      pco_estimate_.dpco_mm_[1] = 0.0;
#endif // WITHOUT_HORIZON_UPDATE
      pcc_->UpdatePCO(pco_estimate_.dpco_mm_);
      // fixしたらpcvのフラグを変える．
      if (pco_estimate_.GetPcoFixed()) pcv_estimate_.SetPcvFixed(false); // estimation_finish_ = true;
      return true;
    }
  }
  else if (!pcv_estimate_.GetPcvFixed())
  {
    if (pcv_estimate_.Update(W, pcc_, elapsed_time))
    {
      pcc_->UpdatePCV(pcv_estimate_.dpcv_vec_mm_);
       // ステップごとに保存できるようにしたい．
      pcc_->PcvLogOutput("_pcv.csv");
      pcc_->PccLogOutput("_pcc.csv");
      if (pcv_estimate_.GetPcvFixed())
      {
#ifdef ITERATION
      // fixしたらpcoのフラグを変えて再度推定を行わせる．
      pco_estimate_.SetPcoFixed(false);
#else
      estimation_finish_ = true;
      pco_estimate_.SetPcoFixed(false); // 初期化しておく
#endif // ITERATION
      }
      return true;
    }
  }
  else std::cout << "ERROR: something wrong" << std::endl;

  return false;
}

