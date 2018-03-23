////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include <cellogram/load_points.h>
#include <cellogram/convex_hull.h>
#include <cellogram/delaunay.h>
#include <cellogram/PolygonUtils.h>
#include <igl/slice.h>
#include <igl/colon.h>
#include <igl/png/readPNG.h>
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

void UIState::initialize() {
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
	viewer.append_mesh();
	viewer.append_mesh();
	viewer.data_list[0].id = hull_id = 0;
	viewer.data_list[1].id = points_id = 1;
	viewer.data_list[2].id = img_id = 2;
}

void UIState::launch() {
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

	// Load points
	load_points(name, state.points);
	state.detected = state.points;

	// Show points and align camera
	points_data().clear();
	points_data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
	viewer.core.align_camera_center(state.points);
	double extent = (state.points.colwise().maxCoeff() - state.points.colwise().minCoeff()).maxCoeff();
	points_data().point_size = float(0.008 * extent);

	// Compute and show convex hull + triangulation
	compute_hull();
	compute_triangulation();

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
	// Compute polygon of the convex hull
	Eigen::MatrixXd P;
	Eigen::VectorXi I, J;
	convex_hull(state.points, state.boundary);
	int dims = (int) state.points.cols();
	I = state.boundary;
	J = Eigen::VectorXi::LinSpaced(dims, 0, dims-1);
	igl::slice(state.points, I, J, P);

	// Offset by epsilon
	// offset_polygon(P, P, 1);

	hull_data().clear();

	// Draw edges
	int n = (int) P.rows();
	Eigen::MatrixXd P2;
	P2.resizeLike(P);
	P2.topRows(n-1) = P.bottomRows(n-1);
	P2.row(n-1) = P.row(0);
	hull_data().add_edges(P, P2, Eigen::RowVector3d(0, 0, 0));

	// Draw filled polygon
	triangulate_convex_polygon(P, state.hull_vertices, state.hull_faces);
	hull_data().set_mesh(state.hull_vertices, state.hull_faces);

	// Set viewer options
	hull_data().set_colors(Eigen::RowVector3d(52, 152, 219)/255.0);
	hull_data().show_lines = false;
	hull_data().shininess = 0;
	hull_data().line_width = 2.0;
}

// -----------------------------------------------------------------------------

void UIState::compute_triangulation() {
	delaunay_triangulation(state.points, state.triangles);
	points_data().clear();
	points_data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
	points_data().set_mesh(state.points, state.triangles);
}

void UIState::load_image(std::string fname) {
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R, G, B, A; // Image

	igl::png::readPNG(fname, R, G, B, A);
	int xMax = R.cols();
	int yMax = R.rows();

	// Replace the mesh with a triangulated square
	Eigen::MatrixXd V(4, 3);
	V <<
		0, 0, 0,
		xMax, 0, 0,
		xMax, yMax, 0,
		0, xMax, 0;
	Eigen::MatrixXi F(2, 3);
	F <<
		0, 1, 2,
		2, 3, 0;
	Eigen::MatrixXd UV(4, 2);
	UV <<
		0, 1,
		1, 1,
		1, 0,
		0, 0;

	img_data().set_mesh(V, F);
	img_data().set_uv(UV);
	img_data().show_texture = true;

	// Use the image as a texture
	img_data().set_texture(R, R, R);

	// Turn of texture of other meshes
	hull_data().show_texture = false;
	points_data().show_texture = false;
}

// -----------------------------------------------------------------------------

} // namespace cellogram
