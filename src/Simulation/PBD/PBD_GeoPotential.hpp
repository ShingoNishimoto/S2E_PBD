#ifndef PBD_GEO_POTENTIAL_H_
#define PBD_GEO_POTENTIAL_H_

#include <vector>
#include <utility>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include <GeoPotential.h>

// 推定系のために重力ポテンシャルを計算する．Dynamicsとかがないので，
class PBD_GeoPotential
{
public:
  PBD_GeoPotential(const int degree, const LocalEnvironment& local_env, std::string file_path);
  ~PBD_GeoPotential();
  Vector<3> CalcAccelerationECI(const Vector<3> &position_eci);

private:
  int degree_;
  const LocalEnvironment& local_env_;
  GeoPotential* geo_potential_;

};


#endif
