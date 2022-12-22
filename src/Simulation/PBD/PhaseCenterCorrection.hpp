#pragma once

#include "../../Library/VectorTool.hpp"
#include <map>

class PhaseCenterCorrection
{
 public:
  PhaseCenterCorrection(libra::Vector<3> pco, std::vector<double> pcv, const double azi_increment, const double ele_increment);
  PhaseCenterCorrection(libra::Vector<3> pco, const double azi_increment, const double ele_increment, const std::string out_fname_base);
  ~PhaseCenterCorrection();

  const double GetPCV_mm(const double azimuth_deg, const double elevation_deg);
  const double GetPCC_m(const double azimuth_deg, const double elevation_deg);
  inline libra::Vector<3> GetPCO_mm(void) const { return pco_mm_; }
  inline void UpdatePCO(const libra::Vector<3> dpco_mm) { pco_mm_ += dpco_mm; }
  void UpdatePCV(const std::vector<double> dpcv_mm);
  inline void SetPCO(const libra::Vector<3> pco_mm) { pco_mm_ = pco_mm; }
  void PcvLogOutput(std::string out_fname);
  void PccLogOutput(std::string out_fname);

  const double azi_increment_;
  const double ele_increment_;
  std::map<int, int> azimuth_index_;   // azimuth角を入力としたindex取得用辞書
  std::map<int, int> elevation_index_; // elevation角を入力としたindex取得用辞書
  std::string out_fname_base_;

 protected:
  libra::Vector<3> pco_mm_;    // 3次元ベクトル
  std::vector<double> pcv_mm_; // grid pointを保持

  void InitAngleIndexes(void);

};
