////////////////////////////////////////////////////////////////////////////////
#include "vertex.h"
#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

void add_vertex(const Eigen::MatrixXd &V){
	// add vertex at clicked position... possibly after fitting a gaussian locally
}
void delete_vertex(const Eigen::MatrixXd &V) {
	// delete vertex closest to clicked position

}

// -----------------------------------------------------------------------------

} // namespace cellogram
