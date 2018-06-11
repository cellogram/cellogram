////////////////////////////////////////////////////////////////////////////////
#include "MeshUtils.h"
#include <iomanip>
#include <cassert>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

void to_geogram_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, GEO::Mesh &M) {
	M.clear();
	// Setup vertices
	M.vertices.create_vertices((int) V.rows());
	for (int i = 0; i < (int) M.vertices.nb(); ++i) {
		GEO::vec3 &p = M.vertices.point(i);
		p[0] = V(i, 0);
		p[1] = V(i, 1);
		p[2] = V(i, 2);
	}
	// Setup faces
	if (F.cols() == 3) {
		M.facets.create_triangles((int) F.rows());
	} else if (F.cols() == 4) {
		M.facets.create_quads((int) F.rows());
	} else {
		throw std::runtime_error("Mesh faces not supported");
	}
	for (int c = 0; c < (int) M.facets.nb(); ++c) {
		for (int lv = 0; lv < F.cols(); ++lv) {
			M.facets.set_vertex(c, lv, F(c, lv));
		}
	}
    M.facets.connect();
}

// -----------------------------------------------------------------------------

void to_geogram_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXi &T, GEO::Mesh &M) {
	to_geogram_mesh(V, F, M);
	if (T.cols() == 4) {
		M.cells.create_tets((int) T.rows());
	} else {
		throw std::runtime_error("Mesh cells not supported");
	}
	for (int c = 0; c < (int) M.cells.nb(); ++c) {
		for (int lv = 0; lv < T.cols(); ++lv) {
			M.cells.set_vertex(c, lv, T(c, lv));
		}
	}
	M.cells.connect();
}

// -----------------------------------------------------------------------------

void from_geogram_mesh(const GEO::Mesh &M, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::MatrixXi &T) {
	V.resize(M.vertices.nb(), 3);
	for (int i = 0; i < (int) M.vertices.nb(); ++i) {
		GEO::vec3 p = M.vertices.point(i);
		V.row(i) << p[0], p[1], p[2];
	}
	assert(M.facets.are_simplices());
	F.resize(M.facets.nb(), 3);
	for (int c = 0; c < (int) M.facets.nb(); ++c) {
		for (int lv = 0; lv < 3; ++lv) {
			F(c, lv) = M.facets.vertex(c, lv);
		}
	}
	assert(M.cells.are_simplices());
	T.resize(M.cells.nb(), 4);
	for (int c = 0; c < (int) M.cells.nb(); ++c) {
		for (int lv = 0; lv < 4; ++lv) {
			T(c, lv) = M.cells.vertex(c, lv);
		}
	}
}

// -----------------------------------------------------------------------------

} // namespace cellogram
