////////////////////////////////////////////////////////////////////////////////
#include "extrude_mesh.h"
#include <igl/boundary_loop.h>
#include <igl/bfs_orient.h>
#include <igl/triangle/triangulate.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
#include <MeshUtils.hpp>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

void extrude_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double thickness, Eigen::MatrixXd &VT, Eigen::MatrixXi &FT, Eigen::MatrixXi &TT) {
	// Extract boundary loop
	Eigen::VectorXi L;
	igl::boundary_loop(F, L);

	// Extract boundary edge-graph
	int n = L.size();
	Eigen::MatrixXd VC(n, 2);
	Eigen::MatrixXi EC(n, 2);
	for (int i = 0; i < n; ++i) {
		VC.row(i) = V.row(L(i)).head<2>();
		EC.row(i) << i, (i+1)%n;
	}

	// Triangulate outer-face
	Eigen::MatrixXd V2;
	Eigen::MatrixXi F2;
	igl::triangle::triangulate(VC, EC, Eigen::MatrixXd(0, 2), "QqY", V2, F2);

	// Sanity check
	Eigen::VectorXi L2;
	igl::boundary_loop(F2, L2);
	cel_assert(L2.size() == n);
	cel_assert(L2.maxCoeff() + 1 == L2.size()); // boundary vertices should be the first indices

	// Stitch both sides
	Eigen::MatrixXd VS(V.rows() + V2.rows(), 3);
	VS.topLeftCorner(V.rows(), 2) = V.leftCols<2>();
	VS.topRightCorner(V.rows(), 1).setZero();
	VS.bottomLeftCorner(V2.rows(), 2) = V2;
	VS.bottomRightCorner(V2.rows(), 1).setConstant(thickness);

	Eigen::MatrixXi FS(F.rows() + F2.rows() + 2 * L.size(), 3);
	FS.topRows(F.rows()) = F;
	FS.middleRows(F.rows(), F2.rows()) = F2.array() + V.rows();
	for (int vo = V.rows(), fo = F.rows() + F2.rows(), i = 0; i < n; ++i) {
		cel_assert(L2(i) == i);
		FS.row(fo + 2*i + 0) << L(i), L((i+1)%n), L2(i)+vo;
		FS.row(fo + 2*i + 1) << L2(i)+vo, L2((i+1)%n)+vo, L((i+1)%n);
	}

	// Orient faces coherently
	Eigen::VectorXi C;
	igl::bfs_orient(FS, FS, C);

	// Mesh volume
	igl::copyleft::tetgen::tetrahedralize(VS, FS, "Qq1.414", VT, TT, FT);
	poly_fem::orient_closed_surface(VT, FT);
}

// -----------------------------------------------------------------------------

} // namespace cellogram

