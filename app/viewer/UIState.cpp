////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"

#include <cellogram/PolygonUtils.h>
#include <igl/colon.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/colormap.h>
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

	//// KEY and MOUSE handling example
	//viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, int key, int mod)->bool
	//{
	//	std::cout << (char)key << " "<< mod << std::endl;
	//	return false;
	//};

	//viewer.callback_key_up = [&](igl::opengl::glfw::Viewer& viewer, int key, int mod)->bool
	//{
	//	std::cout << (char)key << " " << mod << std::endl;
	//	return false;
	//};


	viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
	{
		if (!select_region && !add_vertex && !delete_vertex && !make_vertex_good)
			return false;

		int fid;
		Eigen::Vector3f bc;
		bool found = false;
		double x = viewer.current_mouse_x;
		double y = viewer.core.viewport(3) - viewer.current_mouse_y;
		Eigen::MatrixXd V = t * state.mesh.points + (1 - t) * state.mesh.detected;

		if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
			viewer.core.proj, viewer.core.viewport, V, state.mesh.triangles, fid, bc))
		{
			//std::cout << fid << " " << bc.transpose() << std::endl;
			found = true;
		}

		if (!found)
			return false;


		if (select_region)
		{
			select_region = false;
			
			// set selected_region
			for (int i = 0; i < state.regions.size(); i++)
			{
				for (int  j = 0; j < state.regions[i].region_triangles.rows(); j++)
				{
					if (fid == state.regions[i].region_triangles(j))
					{
						selected_region = i;
						viewer_control();
						current_region_status = Region::pretty_status(state.regions[i].status);
						std::cout << current_region_status << std::endl;
						return true;
					}
				}
			}
		}
		else if(add_vertex)
		{ 
			add_vertex = false;
			double xNew = 0, yNew = 0, zNew = 0;
			for (int i = 0; i < 3; i++)
			{
				xNew += bc(i) * V(state.mesh.triangles(fid, i), 0);
				yNew += bc(i) * V(state.mesh.triangles(fid, i), 1);
				zNew += 0;
			}
			state.add_vertex(Eigen::Vector3d(xNew,yNew,zNew));
			state.reset_state();
			reset_viewer();
			viewer_control();
		}
		else if (delete_vertex)
		{
			delete_vertex = false;
			// find maximum barycenter coordinates and get vertex id
			int vid;
			int maxBC = 0;
			double maxBCVal = bc(maxBC);
			for (int i = 1; i < 3; i++)
			{
				if (maxBCVal < bc(i)) {
					maxBC = i;
					maxBCVal = bc(i);
				}
			}
			//std::cout << "\n max\n" <<  maxBC;
			vid = state.mesh.triangles(fid, maxBC);
			state.delete_vertex(vid);
			//state.reset_state();
			reset_viewer();
			viewer_control();
		}
		else if (make_vertex_good)
		{
			make_vertex_good = false;
			// find maximum barycenter coordinates and get vertex id
			int vid;
			int maxBC = 0;
			double maxBCVal = bc(maxBC);
			for (int i = 1; i < 3; i++)
			{
				if (maxBCVal < bc(i)) {
					maxBC = i;
					maxBCVal = bc(i);
				}
			}
			//std::cout << "\n max\n" <<  maxBC;
			vid = state.mesh.triangles(fid, maxBC);
			state.fixed_as_good.push_back(vid);
		}
		return false;
	};

	viewer.plugins.push_back(this);

	/*state.load("C:\\Users\\Tobias\\Documents\\cellogram\\data\\small2.xyz");
	state.compute_hull();
	state.compute_triangulation();
	state.relax_with_lloyd();
	state.detect_bad_regions();
	state.fix_regions();
	state.resolve_regions();*/

	// Setup viewer parameters
	viewer.resize(1500, 1500);
	viewer.core.background_color.setOnes();
	// viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);
	viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);
	viewer.core.orthographic = true;
	viewer.core.is_animating = true;
	viewer.core.is_animating = true;

	// Setup viewer data
	viewer.append_mesh();
	viewer.append_mesh();
	viewer.append_mesh();
	viewer.append_mesh();
	viewer.data_list[0].id = hull_id = 0;
	viewer.data_list[1].id = points_id = 1;
	viewer.data_list[2].id = image_id = 2;
	viewer.data_list[3].id = bad_region_id = 3;
	viewer.data_list[4].id = matching_id = 4;
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
	if (!state.load(name)) { return false; }
	selected_region = -1;
	current_region_status = "";
	//img.resize(0, 0);
	// reset flags

	mesh_color.resize(1, 3);
	mesh_color.row(0) = Eigen::RowVector3d(255, 255, 120) / 255.0;
	reset_viewer();

	// Show points and align camera
	points_data().clear();
	points_data().set_points(state.mesh.points, Eigen::RowVector3d(1, 0, 0));
	viewer.core.align_camera_center(state.mesh.points);
	double extent = (state.mesh.points.colwise().maxCoeff() - state.mesh.points.colwise().minCoeff()).maxCoeff();
	points_data().point_size = float(0.008 * extent);
	fix_color(points_data());

	// Compute and show convex hull + triangulation
	compute_hull();
	//clean_hull();
	compute_triangulation();
	viewer_control();
	return true;
}

void UIState::detect_vertices() {
	state.detect_vertices();
	selected_region = -1;
	current_region_status = "";
	//img.resize(0, 0);
	// reset flags

	mesh_color.resize(1, 3);
	mesh_color.row(0) = Eigen::RowVector3d(255, 255, 120) / 255.0;
	reset_viewer();

	// Show points and align camera
	points_data().clear();
	points_data().set_points(state.mesh.points, Eigen::RowVector3d(1, 0, 0));
	viewer.core.align_camera_center(state.mesh.points);
	double extent = (state.mesh.points.colwise().maxCoeff() - state.mesh.points.colwise().minCoeff()).maxCoeff();
	points_data().point_size = float(0.008 * extent);
	fix_color(points_data());

	// Compute and show convex hull + triangulation
	compute_hull();
	//clean_hull();
	compute_triangulation();
	viewer_control();
}

bool UIState::load_param(std::string name)
{
	return state.load_param(name);
}

bool UIState::save(std::string name) {
	return state.save(name);
}

bool UIState::mouse_scroll(float delta_y) {
	viewer.scroll_position += delta_y;

	// Only zoom if there's actually a change
	if (delta_y != 0)
	{
		float mult = (1.0 + ((delta_y>0) ? 1. : -1.)*0.1);
		const float min_zoom = 0.1f;
		viewer.core.camera_zoom = (viewer.core.camera_zoom * mult > min_zoom ? viewer.core.camera_zoom * mult : min_zoom);
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////

void UIState::compute_hull() {
	// State
	state.compute_hull();
	
	// UI
	show_hull = true;
}

void UIState::clean_hull() {
	// State
	state.clean_hull();
}

// -----------------------------------------------------------------------------

void UIState::compute_triangulation() {
	// State
	//state.compute_triangulation();

	// UI
	show_points = true;
}

void UIState::reset_viewer()
{
	// Display flags
	t = 0;
	vertex_color = Eigen::RowVector4f(1,0,0,0);
	mesh_color.resize(0, 0);
	show_hull = false;
	show_points = false;
	show_mesh_fill = true;
	show_image = false;
	show_matching = false;
	show_bad_regions = false;
}

//TODO refactor when more clear
void UIState::load_image(std::string fname) {
	state.load_image(fname);

	show_mesh_fill = false;
	show_image = true;

}

void UIState::display_image()
{
	int xMax = state.img.cols();
	int yMax = state.img.rows();

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
		0, 0,
		1, 0,
		1, 1,
		0, 1;

	image_data().set_mesh(V, F);
	image_data().set_uv(UV);
	image_data().show_faces = true;
	image_data().show_lines = false;
	image_data().show_texture = true;
	// Use the image as a texture
	image_data().set_texture(state.img, state.img, state.img);
	image_data().set_colors(Eigen::RowVector3d(1, 1, 1));
}

void UIState::viewer_control()
{
	hull_data().clear();
	points_data().clear();
	image_data().clear();
	matching_data().clear();
	bad_region_data().clear();

	// Read all the UIState flags and update display
	// hull
	if (show_hull) {
		// Draw edges
		int n = (int)state.hull_polygon.rows();

		Eigen::MatrixXd P2; P2.resizeLike(state.hull_polygon);
		P2.topRows(n - 1) = state.hull_polygon.bottomRows(n - 1);
		P2.row(n - 1) = state.hull_polygon.row(0);
		hull_data().add_edges(state.hull_polygon, P2, Eigen::RowVector3d(0, 0, 0));

		hull_data().set_mesh(state.hull_vertices, state.hull_faces);
		hull_data().set_colors(Eigen::RowVector3d(52, 152, 219) / 255.0);
		hull_data().show_faces = false;
		hull_data().show_lines = false;
		hull_data().shininess = 0;
		hull_data().line_width = 3.0;
	}

	// points
	Eigen::MatrixXd V = t * state.mesh.points + (1 - t) * state.mesh.detected;
	points_data().set_mesh(V, state.mesh.triangles);
	/*if (image_loaded)
	{
		Eigen::MatrixXd UV(state.detected.rows(), 2);
		UV.col(0) = state.detected.col(0) / img.rows();
		UV.col(1) = 1 - state.detected.col(1).array() / img.cols();

		points_data().show_texture = true;
		points_data().set_uv(UV);
		points_data().set_texture(img, img, img);
		points_data().set_colors(Eigen::RowVector3d(1, 1, 1));
	}*/
	if (show_points)
	{
		if (state.regions.empty())
			points_data().set_points(V, vertex_color.cast<double>().head<3>());
		else
		{
			Eigen::MatrixXd C(V.rows(), 3);
			C.col(0).setConstant(vertex_color(0));
			C.col(1).setConstant(vertex_color(1));
			C.col(2).setConstant(vertex_color(2));

			if (color_code)
			{
				//Eigen::VectorXd param(V.rows());
				//for (int i = 0; i < V.rows(); i++)
				//{
				//	param(i) = state.mesh.params(i,selected_param);
				//}
				////igl::parula(param, true, C);
				//igl::ColorMapType cm = igl::ColorMapType::COLOR_MAP_TYPE_INFERNO;
				//igl::colormap(cm,param, true, C);
			}
			else
			{
				for (auto &r : state.regions)
				{

					for (int i = 0; i < r.region_boundary.size(); ++i)
					{
						C.row(r.region_boundary(i)) = Eigen::RowVector3d(0, 0, 1);
					}

					for (int i = 0; i < r.region_interior.size(); ++i)
					{
						C.row(r.region_interior(i)) = Eigen::RowVector3d(0, 1, 0);
					}

				}
			}


			points_data().set_points(V, C);

		}
	}

	// fill
	points_data().show_faces = show_mesh_fill;
	if (mesh_color.size() > 0) {
		if (selected_region >= 0) {
			for (int i = 0; i < state.regions[selected_region].region_interior.size(); ++i)
				mesh_color.row(state.regions[selected_region].region_interior[i]) = Eigen::RowVector3d(0, 1, 0);
		}
	}
	points_data().set_colors(mesh_color);

	// bad regions
	if (show_bad_regions)
	{
		Eigen::MatrixXd bad_P1, bad_P2;
		build_region_edges(V, bad_P1, bad_P2);
		bad_region_data().add_edges(bad_P1, bad_P2, Eigen::RowVector3d(0, 0, 0));
		bad_region_data().line_width = 3.0;
	}


	// image
	if (show_image)
	{
		if (state.img.size() <= 0)
		{
			show_image = false;
		}
		else
		{
			display_image();
		}
	}

	// matching
	if (show_matching) {
		matching_data().clear();
		matching_data().add_edges(state.mesh.points, state.mesh.detected, Eigen::RowVector3d(0, 0, 0));
		matching_data().line_width = 3.0;
	}

	// Fix shininess for all layers
	fix_color(hull_data());
	fix_color(points_data());
	fix_color(image_data());
	fix_color(bad_region_data());
	fix_color(matching_data());
}

void UIState::draw_mesh()
{
	points_data().clear();
	points_data().set_points(state.mesh.points, Eigen::RowVector3d(1, 0, 0));
	points_data().set_mesh(state.mesh.points, state.mesh.triangles);

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
	Eigen::VectorXd regionD(state.mesh.points.rows());
	regionD.setZero();
	for (int i = 0; i < state.regions.size(); ++i)
	{
		/*for (int j = 0; j < state.regions[i].region_boundary.size(); ++j)
		{
			regionD(state.regions[i].region_boundary(j)) = i + 10;
		}*/

		for (int j = 0; j < state.regions[i].region_interior.size(); ++j)
		{
			regionD(state.regions[i].region_interior(j)) = i+1;
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
