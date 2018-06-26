#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Computes the laplace energy of a 2d graph}
///
/// @param[in]  V     { #V x dims input point positions }
/// @param[in]  F     { #F x 3 input triangle indices }
/// @param[out] E     { #1 x dims output laplace energy at vertices }
///
void laplace_energy(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXd &E);

} // namespace cellogram
