////////////////////////////////////////////////////////////////////////////////
#include "load_points.h"
#include "MeshUtils.h"
#include <geogram/mesh/mesh_io.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

void load_points(const std::string &filename, Eigen::MatrixXd &V) {
	GEO::Mesh M;
	if(!GEO::mesh_load(filename, M)) {
		return;
	}
	Eigen::MatrixXi F, T;
	from_geogram_mesh(M, V, F, T);
}

} // namespace cellogram
