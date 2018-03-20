////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

std::vector<int> convex_hull(const std::vector<GEO::vec2> &points);

void triangulate_hull(std::vector<GEO::vec2> &hull, GEO::Mesh &M);

///
/// Compute the convex hull of a set of 2D points. Only their XY coordinate is used.
///
/// @param[in]  V     { #V x (2|3) input point positions }
/// @param[out] L     { #L list of indices into V of points that form the convex hull of V. }
///
void convex_hull(const Eigen::MatrixXd &V, Eigen::VectorXi &L);

///
/// Triangulate a convex polygon
///
/// @param[in]  P     { #P x (2|3) input point positions along the polygon }
/// @param[out] V     { #V x (2|3) output vertex positions }
/// @param[out] F     { #F x 3 output triangle indices }
///
void triangulate_convex_polygon(const Eigen::MatrixXd &P, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

} // namespace cellogram
