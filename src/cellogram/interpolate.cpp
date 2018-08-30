////////////////////////////////////////////////////////////////////////////////
#include "interpolate.h"
#include "MeshUtils.h"
#include <geogram/mesh/mesh_AABB.h>
#include <vector>
#include <queue>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

namespace {

	// http://gamedev.stackexchange.com/a/23745
	void barycentric(GEO::vec3 p, GEO::vec3 a, GEO::vec3 b, GEO::vec3 c, double &u, double &v, double &w) {
		GEO::vec3 v0 = b - a, v1 = c - a, v2 = p - a;
		float d00 = dot(v0, v0);
		float d01 = dot(v0, v1);
		float d11 = dot(v1, v1);
		float d20 = dot(v2, v0);
		float d21 = dot(v2, v1);
		float denom = d00 * d11 - d01 * d01;
		v = (d11 * d20 - d01 * d21) / denom;
		w = (d00 * d21 - d01 * d20) / denom;
		u = 1.0f - v - w;
	}

} // anonymous namespace

// -----------------------------------------------------------------------------

Eigen::VectorXd interpolate_2d(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &S, const Eigen::MatrixXd &Q) {
	GEO::Mesh M;
	to_geogram_mesh(V, F, M);
	GEO::MeshFacetsAABB aabb_tree(M, false);
	Eigen::VectorXd X(Q.rows());
	X.setZero();
	for (int q = 0; q < Q.rows(); ++q) {
		GEO::vec3 pts(Q(q,0), Q(q,1), Q.cols() < 3 ? 0 : Q(q,2));
		GEO::vec3 proj;
		double sq_dist;
		GEO::index_t f = aabb_tree.nearest_facet(pts, proj, sq_dist);
		GEO::index_t ia = M.facets.vertex(f, 0); GEO::vec3 a = M.vertices.point(ia);
		GEO::index_t ib = M.facets.vertex(f, 1); GEO::vec3 b = M.vertices.point(ib);
		GEO::index_t ic = M.facets.vertex(f, 2); GEO::vec3 c = M.vertices.point(ic);
		double u, v, w;
		barycentric(proj, a, b, c, u, v, w);
		u = std::max(0.0, std::min(1.0, u));
		v = std::max(0.0, std::min(1.0, v));
		w = std::max(0.0, std::min(1.0, w));
		X(q) = u*S(ia) + v*S(ib) + w*S(ic);
	}
	return X;
}

} // namespace cellogram
