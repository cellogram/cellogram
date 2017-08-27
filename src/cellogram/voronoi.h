////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

void lloyd_relaxation(std::vector<GEO::vec2> &points, const std::vector<bool> &fixed, int num_iter,
	GEO::Mesh *domain = nullptr);

} // namespace cellogram
