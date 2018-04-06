////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"

#include <cellogram/PolygonUtils.h>
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
	if (!state.load(name)){ return false; }

	// Show points and align camera
	points_data().clear();
	points_data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
	viewer.core.align_camera_center(state.points);
	double extent = (state.points.colwise().maxCoeff() - state.points.colwise().minCoeff()).maxCoeff();
	points_data().point_size = float(0.008 * extent);
	fix_color(points_data());

	// Compute and show convex hull + triangulation
	compute_hull();
	compute_triangulation();

	return true;
}

bool UIState::save(std::string name) {
	return state.save(name);
}

////////////////////////////////////////////////////////////////////////////////

void UIState::compute_hull() {
	state.compute_hull();
	


	hull_data().clear();

	// Draw edges
	int n = (int) state.hull_polygon.rows();

	Eigen::MatrixXd P2; P2.resizeLike(state.hull_polygon);
	P2.topRows(n-1) = state.hull_polygon.bottomRows(n-1);
	P2.row(n-1) = state.hull_polygon.row(0);
	hull_data().add_edges(state.hull_polygon, P2, Eigen::RowVector3d(0, 0, 0));

	hull_data().set_mesh(state.hull_vertices, state.hull_faces);

	// Set viewer options
	hull_data().set_colors(Eigen::RowVector3d(52, 152, 219)/255.0);
	hull_data().show_faces = false;
	hull_data().show_lines = false;
	hull_data().shininess = 0;
	hull_data().line_width = 2.0;
}

// -----------------------------------------------------------------------------

void UIState::compute_triangulation() {
	state.compute_triangulation();

	draw_mesh();
}

//TODO refactor when more clear
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
	img_data().show_faces = true;
	img_data().show_lines = false;
	// Use the image as a texture
	img_data().set_texture(R,G,B);
}

void UIState::draw_mesh()
{
	points_data().clear();
	points_data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
	points_data().set_mesh(state.points, state.triangles);

	fix_color(points_data());
}

void UIState::fix_color(igl::opengl::ViewerData &data) {
	data.F_material_specular.setZero();
	data.V_material_specular.setZero();
	data.dirty |= igl::opengl::MeshGL::DIRTY_DIFFUSE;

	data.V_material_ambient *= 2;
	data.F_material_ambient *= 2;
}


Eigen::VectorXd UIState::create_region_label()
{
	Eigen::VectorXd regionD(state.points.rows());
	regionD.setZero();
	for (int i = 0; i < state.regions.size(); ++i)
	{
		/*for (int j = 0; j < state.regions[i].region_boundary.size(); ++j)
		{
			regionD(state.regions[i].region_boundary(j)) = i + 10;
		}*/

		for (int j = 0; j < state.regions[i].region_interior.size(); ++j)
		{
			regionD(state.regions[i].region_interior(j)) = i + 10;
		}
	}

	return regionD;
}

void UIState::build_region_edges(const Eigen::MatrixXd &pts, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2)
{
	bad_P1.resize(pts.rows(), 3);
	bad_P2.resize(pts.rows(), 3);

	int index = 0;
	Eigen::MatrixXd local_bad_P1, local_bad_P2;
	for (auto &r : state.regions)
	{
		r.compute_edges(pts, local_bad_P1, local_bad_P2);
		bad_P1.block(index, 0, local_bad_P1.rows(), local_bad_P1.cols()) = local_bad_P1;
		bad_P2.block(index, 0, local_bad_P2.rows(), local_bad_P2.cols()) = local_bad_P2;

		index += local_bad_P2.rows();

	}
	bad_P1.conservativeResize(index, 3);
	bad_P2.conservativeResize(index, 3);
}
// -----------------------------------------------------------------------------

} // namespace cellogram
