#ifndef __vector_tool_H__
#define __vector_tool_H__

#include <Dense> // Eigen
#include <vector>
#include "Vector.hpp"

inline libra::Vector<3> ConvEigenVecToLibraVec(Eigen::Vector3d eigen_vec)
{
  libra::Vector<3> libra_vec(0);
  for (uint8_t i = 0; i < 3; i++) libra_vec[i] = eigen_vec(i);

  return libra_vec;
}

inline Eigen::Vector3d ConvLibraVecToEigenVec(libra::Vector<3> libra_vec)
{
  Eigen::Vector3d eigen_vec = Eigen::Vector3d::Zero();
  for (uint8_t i = 0; i < 3; i++) eigen_vec(i) = libra_vec[i];

  return eigen_vec;
}

#endif __vector_tool_H__
