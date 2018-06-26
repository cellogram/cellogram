#pragma once

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

///
/// Compute a loose convex hull around the given set of points. This function
/// computes an initial convex hull and Delaunay triangulation of the 2D point
/// cloud, then prunes off boundary triangles with a border edge that is more
/// than `edge_length_ratio` times the median edge-length of the initial
/// triangulation.
///
/// @param[in]  V                  { #V x (2|3) input point positions }
/// @param[out] L                  { #L list of indices into V of points along
///                                the relaxed convex hull. }
/// @param[in]  edge_length_ratio  { Edge length ratio used for the pruning }
///
void loose_convex_hull(const Eigen::MatrixXd &V, Eigen::VectorXi &L, double edge_length_ratio = 3.0);

///
/// Triangulate a (potentially non-convex) polygon
///
/// @param[in]  P     { #P x (2|3) input point positions along the polygon }
/// @param[out] V     { #V x (2|3) output vertex positions }
/// @param[out] F     { #F x 3 output triangle indices }
///
void triangulate_polygon(const Eigen::MatrixXd &P, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
//void triangulate_polygon(const Eigen::MatrixXd &P, std::vector<std::vector<int>> &regions, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

} // namespace cellogram
