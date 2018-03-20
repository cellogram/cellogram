////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <geogram/basic/geometry.h>
#include <geogram/mesh/mesh.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

std::vector<int> convex_hull(const std::vector<GEO::vec2> &points);

void triangulate_hull(std::vector<GEO::vec2> &hull, GEO::Mesh &M);

} // namespace cellogram
