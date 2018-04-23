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

	//put here all helpers functions

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////



void point_source_detection(const Eigen::MatrixXd &img, const double sigma, Eigen::MatrixXd &V)
{
	assert(img.minCoeff()>=0);
	assert(img.maxCoeff()<=1);

	//TODO
}

} // namespace cellogram
