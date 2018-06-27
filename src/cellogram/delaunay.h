#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <cellogram/common.h>
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// Computes a Delaunay tiangulation of a point cloud in 2d
///
/// @param[in]  V     { #V x dims input point positions }
/// @param[out] F     { #F x 3 output triangle indices }
///
void delaunay_triangulation(const Eigen::MatrixXd &V, Eigen::MatrixXi &F);

///
/// Delaunay tiangulation with edge constraints
///
/// @param[in]  V     { #V x dims input point positions }
/// @param[in]  L     { closed polygon forming the boundary of the domain }
/// @param[out] F     { #F x 3 output triangle indices }
///
void constrained_delaunay_triangulation(const Eigen::MatrixXd &V, const Eigen::VectorXi &L, Eigen::MatrixXi &F);

} // namespace cellogram
