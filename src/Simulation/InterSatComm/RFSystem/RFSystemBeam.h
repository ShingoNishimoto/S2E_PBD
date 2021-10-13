#pragma once
#include <string>
#include "Vector.hpp"

// 後で名前修正した方が良さそう
class RFSystemBeam
{
public:
  RFSystemBeam();
  ~RFSystemBeam();

  // こいつらどうする
  //Setter
  inline void SetFreqShift(libra::Vector<2> freq_shift_hz) { freq_shift_hz_ = freq_shift_hz; }
  //Getter
  inline const libra::Vector<2> GetFreqShift(void) const { return freq_shift_hz_; }

private:
  libra::Vector<2> freq_shift_hz_{0.0};
};
