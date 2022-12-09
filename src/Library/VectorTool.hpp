#ifndef __vector_tool_H__
#define __vector_tool_H__

#include <Dense> // Eigen
#include <vector>
#include "Vector.hpp"

template <size_t N>
inline libra::Vector<N> ConvEigenVecToLibraVec(Eigen::VectorXd eigen_vec)
{
  libra::Vector<N> libra_vec(0);
  for (uint8_t i = 0; i < N; i++) libra_vec[i] = eigen_vec(i);

  return libra_vec;
}

template <size_t N>
inline Eigen::VectorXd ConvLibraVecToEigenVec(libra::Vector<N> libra_vec)
{
  // もはや関数にする意味はないかも
  Eigen::VectorXd eigen_vec = Eigen::Map<Eigen::VectorXd>(&libra_vec[0], N);

  return eigen_vec;
}

inline std::vector<double> ConvEigenVecToStdVec(Eigen::VectorXd eigen_vec)
{
  const int len = eigen_vec.rows();
  std::vector<double> std_vec(len); // 要素数指定
  for (uint8_t i = 0; i < len; i++) std_vec.at(i) = eigen_vec(i);

  return std_vec;
}

inline Eigen::VectorXd ConvStdVecToEigenVec(std::vector<double> std_vec)
{
  Eigen::VectorXd eigen_vec = Eigen::Map<Eigen::VectorXd>(&std_vec.at(0), std_vec.size());

  return eigen_vec;
}

// なぜかconst vectorにするとビルド通らない．
template <typename T>
inline const int GetIndexOfStdVector(std::vector<T> target_vector, const T target_value)
{
  std::vector<T>::iterator itr = std::find(target_vector.begin(), target_vector.end(), target_value);
  if (itr == target_vector.end()) // なんかboolを探索するとfalseになる
  {
    std::cout << "not found" << target_value << std::endl;
    // abort();
    return -1;
  }
  const int index = std::distance(target_vector.begin(), itr);
  return index;
};

#endif __vector_tool_H__
