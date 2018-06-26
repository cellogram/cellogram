////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// Extrude a 2D mesh into a volumetric tet-mesh. Assumes the input mesh has no holes.
///
/// @param[in]  V          { #V x 2 input mesh vertices }
/// @param[in]  F          { #F x 3 input mesh triangles }
/// @param[in]  thickness  { Target extruded thickness }
/// @param[out] VT         { #VT x 3 output mesh vertices }
/// @param[out] FT         { #FT x 3 output mesh triangles }
/// @param[out] TT         { #TT x 4 output mesh tetrahedra }
///
void extrude_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double thickness, Eigen::MatrixXd &VT, Eigen::MatrixXi &FT, Eigen::MatrixXi &TT);

} // namespace cellogram
