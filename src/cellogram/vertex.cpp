////////////////////////////////////////////////////////////////////////////////
#include "delaunay.h"
#include <geogram/mesh/mesh.h>
#include <geogram/delaunay/delaunay.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

void add_vertex(const Eigen::MatrixXd &V){
	// add vertex at clicked position... possibly after fitting a gaussian locally
}
void delete_vertex(const Eigen::MatrixXd &V) {
	// delete vertex closest to clicked position
}

// -----------------------------------------------------------------------------

} // namespace cellogram
