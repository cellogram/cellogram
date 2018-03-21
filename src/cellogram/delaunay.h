////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// @brief      { Computes a Delaunay tiangulation of a point cloud in 2d }
///
/// @param[in]  V     { #V x dims input point positions }
/// @param[out] F     { #F x 3 output triangle indices }
///
void delaunay_triangulation(const Eigen::MatrixXd &V, Eigen::MatrixXi &F);

} // namespace cellogram
