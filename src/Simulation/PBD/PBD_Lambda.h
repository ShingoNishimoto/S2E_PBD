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
  PBD_Lambda(Eigen::MatrixXd P_a, Eigen::VectorXd a_est, const int n_cands); // éQè∆ìnÇµÇ…Ç∑Ç◊Ç´Ç≈ÇÕ?
  ~PBD_Lambda();

  void Solve(void);
 
private:
  int n;
  int n_cands;
  Eigen::MatrixXd Z; // Z transformation
  Eigen::MatrixXd L; // Lower triangular matrix
  Eigen::VectorXd a_hat; // estimated ambiguity
  Eigen::VectorXd a; // fixed ambiguity
  Eigen::MatrixXd z; // fixed z transformed ambigiuty
  Eigen::VectorXd sq_norm;
  Eigen::VectorXd D; // Diagnal of LtDL of P_a

  void IntegerGaussTransformation(const int i, const int j);
  void Permute(const int k, const double delta);
  void Reduction(void);
  void IntegerLeastSquares(const int n_cands, Eigen::VectorXd z_hat);
  Eigen::Index Argmax(const Eigen::VectorXd& x);
  Eigen::Index Argmin(const Eigen::VectorXd& x);
  std::vector<Eigen::Index> Argsort(Eigen::VectorXd x);
  void RemoveRow(Eigen::VectorXd& matrix, unsigned int row);
    // std::random_device seed_gen;
  std::mt19937 mt;
};


#endif
