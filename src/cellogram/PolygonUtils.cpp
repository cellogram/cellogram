////////////////////////////////////////////////////////////////////////////////
#include "PolygonUtils.h"
#include <clipper/clipper.hpp>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------


	typedef std::complex<double> P;

	double inline det(const P &u, const P &v) {
		return imag(conj(u) * v);
	}
	// Return true iff [a,b] intersects [c,d], and store the intersection in ans
	bool intersect_segment(const P &a, const P &b, const P &c, const P &d, P &ans) {
		const double eps = 1e-10; // small epsilon for numerical precision
		double x = det(c - a, d - c);
		double y = det(b - a, a - c);
		double z = det(b - a, d - c);
		// ab and cd are parallel ||
		if (std::abs(z) < eps || x * z < 0 || x * z > z*z || y * z < 0 || y * z > z*z) return false;
		ans = c + (d - c) * y / z;
		return true;
	}

	bool is_inside(const Eigen::MatrixXd &poly, double &x, double &y) {
		P outside(-1000, -1000);
		P q(x, y);
		size_t n = poly.rows();
		bool tmp, ans = false;
		for (size_t i = 0; i < poly.rows(); ++i) {
			P m; // Coordinates of intersection point
			P p0(poly(i, 0), poly(i, 1));
			P p1(poly((i + 1) % n, 0), poly((i + 1) % n, 1));
			tmp = intersect_segment(q, outside, p0, p1, m);
			ans = (ans != tmp);
		}
		return ans;
	}


	// -----------------------------------------------------------------------------



namespace {

constexpr int FACTOR = 1<<23;

namespace Point {
	ClipperLib::IntPoint toClipper(const Eigen::RowVector2d& p) {
		ClipperLib::IntPoint r;
		r.X = (ClipperLib::cInt)std::round(p.x() * FACTOR);
		r.Y = (ClipperLib::cInt)std::round(p.y() * FACTOR);
		return r;
	}

	Eigen::RowVector2d fromClipper(const ClipperLib::IntPoint& p) {
		return Eigen::RowVector2d(p.X, p.Y) / FACTOR;
	}
}

ClipperLib::Path toClipper(const Eigen::MatrixXd &V) {
	ClipperLib::Path path(V.rows());
	for (size_t i = 0; i < path.size(); ++i) {
		path[i] = Point::toClipper(V.row(i));
	}
	return path;
}

Eigen::MatrixXd fromClipper(const ClipperLib::Path &path) {
	Eigen::MatrixXd V(path.size(), 2);
	for (size_t i = 0; i < path.size(); ++i) {
		V.row(i) = Point::fromClipper(path[i]);
	}
	return V;
}



} // anonymous namespace

// -----------------------------------------------------------------------------

void offset_polygon(const Eigen::MatrixXd &IV, Eigen::MatrixXd &OV, double eps) {
	using namespace ClipperLib;

	// Convert input polygon to integer grid
	ClipperOffset co;
	co.AddPath(toClipper(IV), jtSquare, etClosedPolygon);

	// Compute offset in the integer grid
	Paths solution;
	co.Execute(solution, eps * (double) FACTOR);
	assert(solution.size() == 1);

	// Convert back to double
	int dims = (int) IV.cols();
	assert(dims == 2 | dims == 3);
	OV = fromClipper(solution.front());
	if (OV.cols() != dims) {
		OV.conservativeResize(OV.rows(), dims);
		OV.col(2).setZero();
	}
}

// -----------------------------------------------------------------------------

} // namespace cellogram
