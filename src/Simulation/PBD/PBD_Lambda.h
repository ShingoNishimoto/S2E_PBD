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
  PBD_Lambda(Eigen::MatrixXd Q_a, Eigen::VectorXd a_est); // éQè∆ìnÇµÇ…Ç∑Ç◊Ç´Ç≈ÇÕ?
  ~PBD_Lambda();

  void Solve(void);
 
private:
  int n;
  Eigen::MatrixXd Z; // Z transformation
  Eigen::MatrixXd L; // Lower triangular matrix
  Eigen::VectorXd a; // estimated ambiguity
  Eigen::VectorXd D; // Diagnal of LtDL of Q_a

  void IntegerGaussTransformation(void);
  void Permute(void);
  // std::random_device seed_gen;
  std::mt19937 mt;
};


#endif
