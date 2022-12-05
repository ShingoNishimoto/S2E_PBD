#ifndef PBD_LAMBDA_H_
#define PBD_LAMBDA_H_

#include <vector>
#include <utility>
#include <random>
#include "Dense" // Eigen

using std::vector;

  // FIXME:データ構造は要検討．この構造にしていいのか？
  struct Ambiguity
  {
    vector<double> N; // [cycle] ambiguity
    vector<int> gnss_sat_id; // 対応するGNSS ID
    vector<bool> is_fixed; //不定性解除に成功したかどうか
  };

class PBD_Lambda
{
public:
  PBD_Lambda(vector<Ambiguity*> N_est_vec, vector<Eigen::MatrixXd>& P_N, const vector<vector<int>> observed_gnss_sat_ids);
  ~PBD_Lambda();

  // まずは全部成功するかそうでないかの0 1にする．
  bool Solve(void);
  // testも追加する．


private:
  int n_;       // 二重差分時の数
  Eigen::MatrixXd Z_; // Z transformation matrix
  Eigen::MatrixXd L_; // Lower triangular matrix
  Eigen::VectorXd N_hat_; // estimated double difference ambiguity
  Eigen::VectorXd N_; // fixed double difference ambiguity
  Eigen::MatrixXd P_ddN_; // double difference ambiguity covariance matrix
  std::pair<int, int> master_d_N_; // masterとなる衛星のindex, d_Nの値を格納．
  Eigen::MatrixXd z_; // fixed z transformed ambiguity
  Eigen::VectorXd sq_norm_;
  Eigen::VectorXd D_; // Diagonal of LtDL of P_a

  vector<Ambiguity*> vec_N_;

  const vector<vector<int>> observed_gnss_sat_ids_; // main, target, commonの情報を格納

  // Nのchと衛星の対応を持っておく必要がある．結局dgpsクラスでやったほうがいい？
  Eigen::MatrixXd InitializeCovariance(vector<Eigen::MatrixXd> P_N);
  void IntegerGaussTransformation(const int i, const int j);
  void Permute(const int k, const double delta);
  void Reduction(void);
  vector<Eigen::VectorXd> IntegerLeastSquares(const int n_cands, Eigen::VectorXd z_hat);
  Eigen::Index Argmax(const Eigen::VectorXd& x);
  Eigen::Index Argmin(const Eigen::VectorXd& x);
  vector<Eigen::Index> Argsort(Eigen::VectorXd x);
  void RemoveRow(Eigen::VectorXd& matrix, unsigned int row);
  bool SuccessRateTest(vector<Eigen::VectorXd> N_candidates);
  bool DiscriminationTest(vector<Eigen::VectorXd> N_candidates);
  void RecoverNoDifference(void);
};


#endif
