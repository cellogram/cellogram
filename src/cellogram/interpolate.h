#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <cellogram/common.h>
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

//
// @brief      { Linear interpolation of a function defined on the vertices of a triangle mesh. }
//
// @param[in]  V     { #V x (2|3) input mesh vertices }
// @param[in]  F     { #F x 3 input mesh facets }
// @param[in]  S     { #V x 1 scalar field to interpolate }
// @param[in]  Q     { #Q x (2|3) query points to interpolate }
//
// @return     { #Q x 1 interpolated data }
//
Eigen::VectorXd interpolate_2d(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &S, const Eigen::MatrixXd &Q);

} // namespace cellogram
