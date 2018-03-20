////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include <cellogram/load_points.h>
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

	viewer.resize(1024, 1024);
	viewer.core.background_color.setOnes();
	viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);
	viewer.core.orthographic = true;
	viewer.core.is_animating = true;
	viewer.core.is_animating = true;

	viewer.launch();
}

bool UIState::load(std::string name) {
	if (name.empty()) { return true; }

	load_points(name, state.points);
	viewer.data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
	viewer.core.align_camera_center(state.points);

	double extent = (state.points.colwise().maxCoeff() - state.points.colwise().minCoeff()).maxCoeff();
	viewer.data().point_size = float(0.008 * extent);

	return false;
}

bool UIState::save(std::string name) {
	if (name.empty()) { return true; }

	if (/*can save*/true) {
		return true;
	}

	return false;
}

// -----------------------------------------------------------------------------

} // namespace cellogram
