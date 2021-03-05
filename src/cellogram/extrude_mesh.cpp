////////////////////////////////////////////////////////////////////////////////
#include "extrude_mesh.h"
#include <polyfem/MeshUtils.hpp>
#include <zebrafish/Logger.hpp>
#include <igl/all_edges.h>
#include <igl/boundary_loop.h>
#include <igl/bfs_orient.h>
#include <igl/orient_outward.h>
#include <igl/write_triangle_mesh.h>
#include <igl/triangle/triangulate.h>
#include <igl/copyleft/tetgen/tetrahedralize.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

using zebrafish::logger;

// -----------------------------------------------------------------------------

void extrude_mesh(const Eigen::MatrixXd &V1, const Eigen::MatrixXi &F1, double thickness, Eigen::MatrixXd &VT, Eigen::MatrixXi &FT, Eigen::MatrixXi &TT) {

	// 1. Create box vertices
	const int n = V1.rows();
	double z0 = (V1.cols() < 3 ? 0 : V1.col(2).mean());
	Eigen::RowVector2d minV = V1.colwise().minCoeff().head<2>();
	Eigen::RowVector2d maxV = V1.colwise().maxCoeff().head<2>();
	Eigen::MatrixXd V2(V1.rows() + 4, 3);
	V2.topLeftCorner(V1.rows(), V1.cols()) = V1;
    V2.col(2).setConstant(z0);
	V2.bottomRows(4) <<
		minV(0), minV(1), z0 + thickness,
		maxV(0), minV(1), z0 + thickness,
		maxV(0), maxV(1), z0 + thickness,
		minV(0), maxV(1), z0 + thickness;
	Eigen::MatrixXi F2 = F1;

	// 2. Project coordinate of boundary points of input surface to lie exactly on the bbox boundary
	Eigen::MatrixXi E;
	igl::all_edges(F1, E);
	std::vector<bool> boundary_vertex(V2.rows(), false);
	double eps = 1e-9 * (maxV - minV).sum();
	for (int e  = 0; e < E.rows(); ++e) {
		for (int lv = 0; lv < 2; ++lv) {
			boundary_vertex[E(e, lv)] = true;
			for (int d  : {0, 1}) {
				double &x = V2(E(e, lv), d);
				if (std::abs(x - minV(d)) < eps) {
					x = minV(d);
				} else if (std::abs(x - maxV(d)) < eps) {
					x = maxV(d);
				}
			}
		}
	}
	for (int i = 0; i < 4; ++i) {
		boundary_vertex.rbegin()[i] = true;
	}

	// 3. Trivial triangulation of boundary vertices
	for (int d = 0; d < 2; ++d) {
		for (double coord : {minV(d), maxV(d)}) {
			std::vector<std::pair<double, int>> top;
			std::vector<std::pair<double, int>> bottom;
			for (int v = 0; v < V2.rows(); ++v) {
				if (boundary_vertex[v] && V2(v, d) == coord) {
					if (v < n) {
						top.emplace_back(V2(v, 1 - d), v);
					} else {
						bottom.emplace_back(V2(v, 1 - d), v);
					}
				}
			}
			cel_assert(bottom.size() == 2);
			std::sort(top.begin(), top.end());
			std::sort(bottom.begin(), bottom.end());
			const int f0 = F2.rows();
			F2.conservativeResize(F2.rows() + top.size(), F2.cols());
			for (size_t i = 0; i < top.size(); ++i) {
				F2.row(f0 + i) << bottom.front().second, top[i].second, (i + 1 == top.size() ? bottom.back().second : top[i + 1].second);
			}
		}
	}

	// 4. Trivial trangulation of bottom side
	const int f0 = F2.rows();
	F2.conservativeResize(F2.rows() + 2, F2.cols());
	F2.row(f0+0) << n+0, n+1, n+2;
	F2.row(f0+1) << n+2, n+0, n+3;

	// 5. Orient facets and tetrahedralize
	Eigen::VectorXi C;
	igl::bfs_orient(F2, F2, C);
	igl::write_triangle_mesh("/Users/ziyizhang/Projects/tmp/debug_before_tetgen.obj", V2, F2);
    std::cerr << "/Users/ziyizhang/Projects/tmp/debug_before_tetgen.obj" << std::endl;
	igl::copyleft::tetgen::tetrahedralize(V2, F2, "Qq1.414", VT, TT, FT);
	polyfem::orient_closed_surface(VT, FT);

	logger().info(" Extrude_mesh: #input_V={}, #input_F={}, thickness={}, #VT={}, #FT={}, #TT={}", n, F1.rows(), thickness, VT.rows(), FT.rows(), TT.rows());
}

void extrude_mesh_old(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, double thickness, Eigen::MatrixXd &VT, Eigen::MatrixXi &FT, Eigen::MatrixXi &TT) {
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
	polyfem::orient_closed_surface(VT, FT);
}

// -----------------------------------------------------------------------------

} // namespace cellogram

