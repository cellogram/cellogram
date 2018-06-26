#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Computes the degree of a graph}
///
/// @param[in]  F		{ #F x 3 input triangle indices }
/// @param[out] degree  { #1 x dims degree of vertices }
///
void vertex_degree(const Eigen::MatrixXi &F, Eigen::VectorXi &degree);

} // namespace cellogram
