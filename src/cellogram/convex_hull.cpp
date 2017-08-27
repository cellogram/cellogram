#include "convex_hull.h"
#include <algorithm>
#include <numeric>
#include <geogram/basic/geometry.h>

////////////////////////////////////////////////////////////////////////////////

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

std::vector<int> cellogram::convex_hull(const std::vector<GEO::vec2> &points) {
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

void cellogram::triangulate_hull(std::vector<GEO::vec2> &hull, GEO::Mesh &M) {
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
