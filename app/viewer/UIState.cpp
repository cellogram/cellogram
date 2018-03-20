////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include <cellogram/load_points.h>
#include <cellogram/convex_hull.h>
#include <igl/slice.h>
#include <igl/colon.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

UIState::UIState()
	: state(State::state())
{ }

UIState &UIState::ui_state() {
	static UIState instance;
	return instance;
}

void UIState::launch() {
	viewer.plugins.push_back(this);

	// Setup viewer parameters
	viewer.resize(1024, 1024);
	viewer.core.background_color.setOnes();
	// viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
	viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);
	viewer.core.orthographic = true;
	viewer.core.is_animating = true;
	viewer.core.is_animating = true;

	// Setup viewer data
	hull_id = viewer.data().id;
	points_id = viewer.append_mesh();

	// Launch viewer
	viewer.launch();
}

////////////////////////////////////////////////////////////////////////////////

igl::opengl::ViewerData & UIState::mesh_by_id(int id) {
	size_t index = viewer.mesh_index(id);
	assert(viewer.data_list[index].id == id);
	return viewer.data_list[index];
}

bool UIState::load(std::string name) {
	if (name.empty()) { return true; }

	load_points(name, state.points);
	points_data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
	viewer.core.align_camera_center(state.points);

	double extent = (state.points.colwise().maxCoeff() - state.points.colwise().minCoeff()).maxCoeff();
	points_data().point_size = float(0.008 * extent);

	compute_hull();

	return false;
}

bool UIState::save(std::string name) {
	if (name.empty()) { return true; }

	if (/*can save*/true) {
		return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////

void UIState::compute_hull() {
	Eigen::VectorXi I, J;
	convex_hull(state.points, I);
	int dims = (int) state.points.cols();
	J = Eigen::VectorXi::LinSpaced(dims, 0, dims-1);
	igl::slice(state.points, I, J, state.hull);

	int n = (int) state.hull.rows();
	Eigen::MatrixXd P;
	P.resizeLike(state.hull);
	P.topRows(n-1) = state.hull.bottomRows(n-1);
	P.row(n-1) = state.hull.row(0);
	hull_data().add_edges(state.hull, P, Eigen::RowVector3d(0, 0, 0));

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	triangulate_convex_polygon(state.hull, V, F);
	hull_data().set_mesh(V, F);
	hull_data().set_colors(Eigen::RowVector3d(52, 152, 219)/255.0);
	hull_data().show_lines = false;
	hull_data().shininess = 0;

	hull_data().line_width = 2.0;
}

// -----------------------------------------------------------------------------

} // namespace cellogram
