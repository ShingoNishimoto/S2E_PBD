#pragma once

#include "VectorTool.hpp"
#include <map>

class PhaseCenterCorrection{
 public:
  PhaseCenterCorrection(libra::Vector<3> pco, std::vector<double> pcv, const double azi_increment, const double ele_increment);

  const libra::Vector<3> GetPCC(const double azimuth_deg, const double elevation_deg);

 protected:
  std::map<int, int> azimuth_index_;   // azimuth角を入力としたindex取得用辞書
  std::map<int, int> elevation_index_; // azimuth角を入力としたindex取得用辞書
  const double azi_increment_;
  const double ele_increment_;
  libra::Vector<3> pco_;
  std::vector<double> pcv_;

};
