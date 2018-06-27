////////////////////////////////////////////////////////////////////////////////
#include "voronoi.h"
#include <cellogram/MeshUtils.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/voronoi/CVT.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

namespace {

void create_bbox_mesh(GEO::vec2 xy_min, GEO::vec2 xy_max, GEO::Mesh &M) {
	M.clear();
	M.vertices.create_vertices(4);
	M.vertices.point(0) = GEO::vec3(xy_min[0], xy_min[1], 0);
	M.vertices.point(1) = GEO::vec3(xy_max[0], xy_min[1], 0);
	M.vertices.point(2) = GEO::vec3(xy_max[0], xy_max[1], 0);
	M.vertices.point(3) = GEO::vec3(xy_min[0], xy_max[1], 0);
	M.facets.create_quad(0, 1, 2, 3);
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

void lloyd_relaxation(
	std::vector<GEO::vec2> &points, const std::vector<bool> &fixed, int num_iter,
	const GEO::Mesh *domain)
{
	using namespace GEO;
	GEO::vec2 xy_min = points[0];
	GEO::vec2 xy_max = points[0];
	for (const GEO::vec2 &p : points) {
		xy_min[0] = std::min(xy_min[0], p[0]);
		xy_max[0] = std::max(xy_max[0], p[0]);
		xy_min[1] = std::min(xy_min[1], p[1]);
		xy_max[1] = std::max(xy_max[1], p[1]);
	}
	GEO::Mesh M;
	if (!domain) {
		create_bbox_mesh(xy_min, xy_max, M);
	} else {
		M.copy(*domain);
	}

	GEO::CentroidalVoronoiTesselation cvt(&M, 2);
	cvt.set_points((int) points.size(), reinterpret_cast<const double *>(&points[0]));
	for (size_t i = 0; i < fixed.size(); ++i) {
		if (fixed[i]) {
			cvt.lock_point((int) i);
		}
	}

	cvt.Lloyd_iterations(num_iter);
	if (num_iter > 10) {
		cvt.Newton_iterations(num_iter);
	}
	std::copy_n(cvt.embedding(0), 2*points.size(), reinterpret_cast<double *>(&points[0]));
}

////////////////////////////////////////////////////////////////////////////////

void lloyd_relaxation(Eigen::MatrixXd &P, const Eigen::VectorXi &fixed, int num_iter,
	const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
	assert(P.cols() == 2 || P.cols() == 3);
	std::vector<GEO::vec2> pts(P.rows());
	std::vector<bool> fixed2(P.rows(), false);
	for (int i = 0; i < P.rows(); ++i) {
		pts[i] = GEO::vec2(P(i, 0), P(i, 1));
	}
	for (int k = 0; k < fixed.size(); ++k) {
		fixed2[fixed[k]] = true;
	}
	GEO::Mesh M;
	to_geogram_mesh(V, F, M);
	lloyd_relaxation(pts, fixed2, num_iter, &M);
	for (int i = 0; i < P.rows(); ++i) {
		P.row(i).head<2>() << pts[i][0], pts[i][1];
	}
}

// -----------------------------------------------------------------------------

} // namespace cellogram
