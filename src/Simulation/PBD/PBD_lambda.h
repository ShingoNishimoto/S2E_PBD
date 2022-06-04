#ifndef PBD_LAMBDA_H_
#define PBD_LAMBDA_H_

#include <vector>
#include <utility>
#include <random>
#include "Dense" // Eigen
#define _USE_MATH_DEFINES
#include <math.h>

class PBD_Lambda
{
public:
  PBD_Lambda(Eigen::MatrixXd Qa, Eigen::VectorXd a_est);
  ~PBD_Lambda();

  void Update();
  void UpdateLambda();
  void CalcIonfreeObservation();
  void UpdateInfoAfterObserved();

private:
  int n;
  Eigen::MatrixXd Z; // Z transformation
  Eigen::MatrixXd L; // Lower triangular matrix

  void IntegerGaussTransformation();

  // std::random_device seed_gen;
  std::mt19937 mt;
};


#endif
