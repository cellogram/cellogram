////////////////////////////////////////////////////////////////////////////////
#include "cellogram/voronoi.h"
#include <geogram/mesh/mesh.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/voronoi/CVT.h>
////////////////////////////////////////////////////////////////////////////////

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

void cellogram::lloyd_relaxation(
	std::vector<GEO::vec2> &points, const std::vector<bool> &fixed, int num_iter,
	GEO::Mesh *domain)
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
		domain = &M;
	}

	GEO::CentroidalVoronoiTesselation cvt(&M, 2);
	cvt.set_points((int) points.size(), reinterpret_cast<const double *>(&points[0]));
	for (size_t i = 0; i < fixed.size(); ++i) {
		if (fixed[i]) {
			cvt.lock_point((int) i);
		}
	}
	cvt.Lloyd_iterations(num_iter);
	//cvt.Newton_iterations(num_iter);
	std::copy_n(cvt.embedding(0), 2*points.size(), reinterpret_cast<double *>(&points[0]));
}
