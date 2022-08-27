#pragma once
#include "OBC.h"
#include "Vector.hpp"

class OBC_Sat0 : public OBC
{
public:
  OBC_Sat0(ClockGenerator* clock_gen);
  ~OBC_Sat0();
  void Initialize();

private:
  void MainRoutine(int count);
};
