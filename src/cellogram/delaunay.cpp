////////////////////////////////////////////////////////////////////////////////
#include "delaunay.h"
#include <igl/triangle/triangulate.h>
#include <igl/write_triangle_mesh.h>
#include <geogram/mesh/mesh.h>
#include <geogram/delaunay/delaunay.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

void delaunay_triangulation(const Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	F.resize(0, 0);
	assert(V.cols() == 2 || V.cols() == 3);
	int n = (int) V.rows();
	const Eigen::MatrixXd P = V.leftCols<2>().transpose();

	// Compute triangulation
	GEO::Delaunay_var delaunay = GEO::Delaunay::create(2, "BDEL2d");
	delaunay->set_vertices(n, P.data());

	// Extract triangles
	F.resize(delaunay->nb_cells(), 3);
	for (int c = 0; c < (int) delaunay->nb_cells(); ++c) {
		for (int lv = 0; lv < 3; ++lv) {
			GEO::signed_index_t v = delaunay->cell_vertex(c, lv);
			F(c, lv) = v;
		}
	}
}

// -----------------------------------------------------------------------------

void constrained_delaunay_triangulation(const Eigen::MatrixXd &V, const Eigen::VectorXi &L, Eigen::MatrixXi &F) {
	if (L.rows() == 0) {
		assert(false);
		delaunay_triangulation(V, F);
		return;
	}
	Eigen::MatrixXd OV, IV = V.leftCols<2>();
	Eigen::MatrixXi E;
	E.resize(L.rows(), 2);
	int m = (int) L.rows();
	E.col(0) = L;
	E.col(1).head(m-1) = L.tail(m-1);
	E.col(1).tail<1>() = L.head<1>();
	igl::triangle::triangulate(IV, E, Eigen::MatrixXd(0, 2), "QpYY", OV, F);
	assert(OV.rows() == V.rows());
}

// -----------------------------------------------------------------------------

} // namespace cellogram
