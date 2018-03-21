////////////////////////////////////////////////////////////////////////////////
#include "PolygonUtils.h"
#include <clipper/clipper.hpp>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

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
