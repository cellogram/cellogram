#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Computes the boundary around a mesh with holes }
///
/// @param[in]  F     { #F x 3 input triangle indices }
/// @param[out] L     { #1 x dims output list of longest boundary }
///
void boundary_loop(const Eigen::MatrixXi &F, Eigen::VectorXi &L);

} // namespace cellogram
