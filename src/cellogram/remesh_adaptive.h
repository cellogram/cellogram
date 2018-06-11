////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

///
/// Remesh a 3d tet-mesh adaptively following to the given scalar field.
///
/// @param[in]  V     { #V x (2|3) input mesh vertices }
/// @param[in]  T     { #T x 4 input mesh tetrahedra }
/// @param[in]  S     { #V x 1 per-vertex scalar field to follow }
/// @param[out] OV    { #OV x (2|3) output mesh vertices }
/// @param[out] OF    { #OF x F output mesh triangles }
/// @param[out] OT    { #OT x 4 output mesh tetrahedra }
///
void remesh_adaptive_3d(const Eigen::MatrixXd &V, const Eigen::MatrixXi &T, const Eigen::VectorXd &S,
	Eigen::MatrixXd &OV, Eigen::MatrixXi &OF, Eigen::MatrixXi &OT);

} // namespace cellogram
