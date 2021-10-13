#pragma once
#include "OBC.h"
#include "Vector.hpp"

class PBD_Components;

class OBC_Sat1 : public OBC
{
public:
  OBC_Sat1(ClockGenerator* clock_gen, PBD_Components& components);
  ~OBC_Sat1();
  void Initialize();

private:
  PBD_Components& components_;
  void MainRoutine(int count);
};
