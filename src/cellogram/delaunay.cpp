////////////////////////////////////////////////////////////////////////////////
#include "delaunay.h"
#include <geogram/mesh/mesh.h>
#include <geogram/delaunay/delaunay.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

namespace {

std::vector<std::pair<int, int>> delaunay_edges(const GEO::Delaunay_var &delaunay) {
	const int dim = delaunay->dimension();

	std::vector<std::pair<int, int>> edges;
	for (int c = 0; c < (int) delaunay->nb_cells(); ++c) {
		for (int lx = 0; lx < dim+1; ++lx) {
			GEO::signed_index_t v1 = delaunay->cell_vertex(c, lx);
			// Triangles in 2D, Tets in 3D -> Each vertex is connected to every
			// other vertex in the cell
			for (int ly = lx+1; ly < dim+1; ++ly) {
				GEO::signed_index_t v2 = delaunay->cell_vertex(c, ly);
				edges.emplace_back(v1, v2);
				edges.emplace_back(v2, v1);
			}
			// tfm::printf("delaunay | e: %s - %s\n", v1, v2);
		}
	}

	// Remove duplicates
	std::sort(edges.begin(), edges.end());
	auto it = std::unique(edges.begin(), edges.end());
	edges.resize(std::distance(edges.begin(), it));

	return edges;
}

} // anonymous namespace

// -----------------------------------------------------------------------------

void delaunay_triangulation(const Eigen::MatrixXd &V, Eigen::MatrixXi &F) {
	assert(V.cols() == 2 | V.cols() == 3);
	int n = (int) V.rows();
	Eigen::MatrixXd P(2, n);
	P = V.leftCols<2>().transpose();

	// Compute triangulation
	GEO::Delaunay_var delaunay = GEO::Delaunay::create(2, "BDEL2d");
	delaunay->set_vertices(n, P.data());

	// Extract triangles
	F.resize(delaunay->nb_cells(), 3);
	for (int c = 0; c < (int) delaunay->nb_cells(); ++c) {
		for (int lv = 0; lv < 3; ++lv) {
			GEO::signed_index_t v = delaunay->cell_vertex(c, lv);
			F(c, lv) = v;
		}
	}
}

// -----------------------------------------------------------------------------

} // namespace cellogram
