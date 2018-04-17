////////////////////////////////////////////////////////////////////////////////
#include "convex_hull.h"
#include "delaunay.h"
#include "navigation.h"
#include <igl/edges.h>
#include <igl/boundary_loop.h>
#include <igl/triangle/cdt.h>
#include <algorithm>
#include <numeric>
#include <stack>
#include <geogram/basic/geometry.h>
#undef IGL_STATIC_LIBRARY
#include <igl/edge_lengths.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

namespace {

struct Compare {
	const std::vector<GEO::vec2> &points;
	int leftmost; // Leftmost point of the poly

	Compare(const std::vector<GEO::vec2> &points_) : points(points_) { }

	bool operator ()(int i1, int i2) {
		if (i2 == leftmost) { return false; }
		if (i1 == leftmost) { return true; }
		GEO::vec2 u = points[i1] - points[leftmost];
		GEO::vec2 v = points[i2] - points[leftmost];
		double d = det(u, v);
		return (d < 0 || (d == 0 && u.length2() < v.length2()));
	}
};

bool inline salientAngle(const GEO::vec2 &a, const GEO::vec2 &b, const GEO::vec2 &c) {
	return (det(b - a, c - a) >= 0);
}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

std::vector<int> convex_hull(const std::vector<GEO::vec2> &points) {
	Compare order(points);
	order.leftmost = 0;
	for(int i = 1; i < (int) points.size(); ++i) {
		if (points[i][0] < points[order.leftmost][0]) {
			order.leftmost = i;
		} else if (points[i][0] == points[order.leftmost][0]
			&& points[i][1] < points[order.leftmost][1])
		{
			order.leftmost = i;
		}
	}
	std::vector<int> index(points.size());
	std::iota(index.begin(), index.end(), 0);
	std::sort(index.begin(), index.end(), order);
	std::vector<int> hull;
	for (int i : index) {
		hull.push_back(i);
		while (hull.size() > 3u && salientAngle(points[hull.end()[-3]],
			points[hull.end()[-2]], points[hull.end()[-1]]))
		{
			hull.end()[-2] = hull.back();
			hull.pop_back(); // Pop inner vertices
		}
	}
	return hull;
}

////////////////////////////////////////////////////////////////////////////////

void triangulate_hull(std::vector<GEO::vec2> &hull, GEO::Mesh &M) {
	auto n = (GEO::index_t) hull.size();
	M.clear();
	M.vertices.create_vertices((int) hull.size());
	for (GEO::index_t i = 0; i < n; ++i) {
		M.vertices.point(i) = GEO::vec3(hull[i][0], hull[i][1], 0);
	}
	M.facets.create_triangles(n - 2);
	for (GEO::index_t i = 1; i + 1 < n; ++i) {
		M.facets.set_vertex(i-1, 0, 0);
		M.facets.set_vertex(i-1, 1, i);
		M.facets.set_vertex(i-1, 2, i+1);
	}
}

////////////////////////////////////////////////////////////////////////////////

void convex_hull(const Eigen::MatrixXd &V, Eigen::VectorXi &P) {
	std::vector<GEO::vec2> pts(V.rows());
	for (int i = 0; i < V.rows(); ++i) {
		pts[i] = GEO::vec2(V(i, 0), V(i, 1));
	}
	auto hull = convex_hull(pts);
	P.resize(hull.size());
	for (int i = 0; i < P.size(); ++i) {
		P(i) = hull[i];
	}
	P.reverseInPlace();
}

////////////////////////////////////////////////////////////////////////////////

void triangulate_convex_polygon(const Eigen::MatrixXd &P, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	int n = (int) P.rows();
	V = P;
	F.resize(n - 2, 3);
	for (int i = 1; i + 1 < n; ++i) {
		F.row(i-1) << 0, i, i+1;
	}
}

////////////////////////////////////////////////////////////////////////////////

void loose_convex_hull(const Eigen::MatrixXd &V, Eigen::VectorXi &L, double edge_length_ratio) {
	// Compute initial Delaunay triangulation
	Eigen::MatrixXi F;
	delaunay_triangulation(V, F);

	// Compute median edge length
	Eigen::MatrixXi E;
	Eigen::VectorXd lengths;
	igl::edges(F, E);
	igl::edge_lengths(V, E, lengths);
	int num_edges = (int) lengths.size();
	std::nth_element(lengths.data(), lengths.data() + num_edges/2, lengths.data() + num_edges);
	double median_edge_length = lengths[num_edges/2];

	// List border edges, prune incident triangle if length is higher than threshold
	NavigationData data(F);
	std::vector<bool> should_remove_face(F.rows(), false); // mark removed triangles
	std::vector<bool> marked_edges(num_edges, false); // mark edges in the pending queue
	std::stack<NavigationIndex> pending_edges;
	for (int f = 0; f < F.rows(); ++f) {
		for (int lv = 0; lv < F.cols(); ++lv) {
			auto index = index_from_face(F, data, f, lv);
			assert(index.edge < num_edges);
			if (switch_face(data, index).face < 0) {
				pending_edges.push(index);
				marked_edges[index.edge] = true;
			}
		}
	}
	while (!pending_edges.empty()) {
		NavigationIndex index = pending_edges.top();
		pending_edges.pop();
		double len = (V.row(index.vertex) - V.row(switch_vertex(data, index).vertex)).norm();
		if (len > edge_length_ratio * median_edge_length) {
			// Remove facet, and visit face-adjacent edges
			should_remove_face[index.face] = true;
			NavigationIndex index2 = index;
			for (int lv = 0; lv < F.cols(); ++lv, index2 = next_around_face(data, index2)) {
				auto index3 = switch_face(data, index2);
				if (index3.face >= 0 && !marked_edges[index2.edge]) {
					pending_edges.push(index3);
					marked_edges[index3.edge] = true;
				}
			}
		}
	}

	// Extract selected faces
	int num_selected = 0;
	for (bool b : should_remove_face) { if (!b) ++num_selected;  }
	Eigen::MatrixXi F2(num_selected, F.cols());
	for (int f = 0, f2 = 0; f < F.rows(); ++f) {
		if (!should_remove_face[f]) {
			F2.row(f2) = F.row(f);
			++f2;
		}
	}

	// Return boundary loop in the reduced triangulation
	igl::boundary_loop(F2, L);
}

////////////////////////////////////////////////////////////////////////////////

void triangulate_polygon(const Eigen::MatrixXd &P, Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	Eigen::MatrixXd PV = P.leftCols<2>();
	Eigen::MatrixXi E, WE;
	Eigen::VectorXi J;
	E.resize(P.rows(), 2);
	int n = (int) P.rows();
	for (int i = 0; i < n; ++i) {
		E.row(i) << i, (i+1)%n;
	}
	igl::triangle::cdt(PV, E, "Q", V, F, WE, J);
	if (P.cols() == 3) {
		V.conservativeResize(V.rows(), 3);
	}
}

//void triangulate_polygon(std::vector<int> &P,  Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
//	Eigen::MatrixXd
//	Eigen::MatrixXd PV = P.leftCols<2>();
//	Eigen::MatrixXi E, WE;
//	Eigen::VectorXi J;
//	E.resize(P.rows(), 2);
//	int n = (int)P.rows();
//	for (int i = 0; i < n; ++i) {
//		E.row(i) << i, (i + 1) % n;
//	}
//	igl::triangle::cdt(PV, E, "V", V, F, WE, J);
//	if (P.cols() == 3) {
//		V.conservativeResize(V.rows(), 3);
//	}
//}
// -----------------------------------------------------------------------------

} // namespace cellogram
