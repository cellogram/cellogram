#pragma once

#include <Eigen/Dense>
#include <vector>

namespace cellogram {

///
/// Compute offset polygon
///
/// @param[in]  IV    { #IV x 2 of vertex positions for the input polygon }
/// @param[out] OV    { #OV x 2 of vertex positions for the offset polygon }
/// @param[in]  eps   { Offset distance }
///
void offset_polygon(const Eigen::MatrixXd &IV, Eigen::MatrixXd &OV, double eps);
bool is_inside(const Eigen::MatrixXd &poly, double &x, double &y);

} // namespace poly_fem
