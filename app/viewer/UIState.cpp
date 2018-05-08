////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"

#include <cellogram/PolygonUtils.h>
#include <cellogram/StringUtils.h>
#include <igl/colon.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/colormap.h>


#include <GLFW/glfw3.h>

#include <sys/types.h> // required for stat.h
#include <sys/stat.h> // no clue why required -- man pages say so

////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace {
		int cellogram_mkdir(const std::string &path)
		{
			int nError;
#if defined(_WIN32)
			std::wstring widestr = std::wstring(path.begin(), path.end());
			nError = _wmkdir(widestr.c_str()); // can be used on Windows
#else
			nError = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // can be used on non-Windows
#endif
			return nError;
		}
	}
// -----------------------------------------------------------------------------

UIState::UIState()
	: state(State::state())
{
	reset_viewer();
}

UIState &UIState::ui_state() {
	static UIState instance;
	return instance;
}

bool UIState::mouse_move(int button, int modifier)
{
	if (dragging_id == -1)
		return false;

	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;

	//reset update position with x and y, store prevx and y to compute delatax, (use current zoom to move point accordingly???)
	int fid;
	Eigen::Vector3f bc;
	igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
		viewer.core.proj, viewer.core.viewport, img_V, img_F, fid, bc);

	double xNew = 0, yNew = 0, zNew = 0;
	for (int i = 0; i < 3; i++)
	{
		//xNew += bc(i) * V(state.mesh.triangles(fid, i), 0);
		xNew += bc(i) * img_V(img_F(fid, i), 0);
		yNew += bc(i) * img_V(img_F(fid, i), 1);
		zNew += 0;
	}

	state.mesh.moved.row(dragging_id) = RowVector3d(xNew, yNew, zNew);
	

	viewer_control();

	return true;
}

bool UIState::mouse_up(int button, int modifier)
{
	if (button != 0 || dragging_id == -1)
		return false;
	//move_vertex = false;
	dragging_id = -1;

	state.reset_state();
	viewer_control();

	return true;
}

bool UIState::mouse_down(int button, int modifier) {
	if (ImGuiMenu::mouse_down(button, modifier)) {
		return true;
	}


	if (button != 0)
		return false;

	//if (!select_region && !add_vertex && !delete_vertex && !make_vertex_good && !make_vertex_bad) {

	//	if (button == 0)
	//	{
	//		this->viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::None;
	//		return true;
	//	}

	//	return false;
	//}

	int fid;
	Eigen::Vector3f bc;
	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;
	Eigen::MatrixXd V = t * state.mesh.moved + (1 - t) * state.mesh.moved;

	if (V.size() <= 0)
		return false;

	if (select_region)
	{
		if (!igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
			viewer.core.proj, viewer.core.viewport, V, state.mesh.triangles, fid, bc)) {
			if (button == 0)
			{
				this->viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::None;
				return true;
			}

			return false;
		}

		select_region = false;

		// set selected_region
		for (int i = 0; i < state.regions.size(); i++)
		{
			for (int j = 0; j < state.regions[i].region_triangles.rows(); j++)
			{
				if (fid == state.regions[i].region_triangles(j))
				{
					selected_region = i;
					viewer_control();
					current_region_status = Region::pretty_status(state.regions[i].status);
					return true;
				}
			}
		}
	}
	else
	{
		if (!igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
			viewer.core.proj, viewer.core.viewport, img_V, img_F, fid, bc))
		{
			if (button == 0)
			{
				this->viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::None;
				return true;
			}

			return false;
		}


		double xNew = 0, yNew = 0, zNew = 0;
		for (int i = 0; i < 3; i++)
		{
			//xNew += bc(i) * V(state.mesh.triangles(fid, i), 0);
			xNew += bc(i) * img_V(img_F(fid, i), 0);
			yNew += bc(i) * img_V(img_F(fid, i), 1);
			zNew += 0;
		}

		Eigen::VectorXd delta;
		delta = (V.col(0).array() - xNew).cwiseAbs2() + (V.col(1).array() - yNew).cwiseAbs2();
		int vid;
		delta.minCoeff(&vid);

		//add_vertex = false;
		if (add_vertex)
		{
			state.add_vertex(Eigen::Vector3d(xNew, yNew, zNew));

			state.reset_state();
			reset_viewer();
			viewer_control();
			return true;
		}
		else if (delete_vertex)
		{
			// find vertex in V closest to xNew,yNew
			state.delete_vertex(vid);

			state.reset_state();
			reset_viewer();
			viewer_control();
			return true;
		}	
		else if (move_vertex)
		{
			dragging_id = vid;
			return true;
		}
		else if (make_vertex_good || make_vertex_bad)
		{
			// find maximum barycenter coordinates and get vertex id
			if (make_vertex_good)
			{
				state.mesh.vertex_status_fixed(vid) = 1;
			}
			else if (make_vertex_bad)
			{
				state.mesh.vertex_status_fixed(vid) = -1;
			}

			make_vertex_good = false;
			make_vertex_bad = false;

			state.detect_bad_regions();
			state.check_regions();
			create_region_label();
			viewer_control();

			return true;
		}
	}

	if (button == 0)
	{
		this->viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::None;
		return true;
	}
	return false;

}

void UIState::initialize() {
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

//bool UIState::load(std::string name) {
//	assert(false);
//	if (!state.load(name)) { return false; }
//	selected_region = -1;
//	current_region_status = "";
//	//img.resize(0, 0);
//	// reset flags
//
//	mesh_color.resize(1, 3);
//	mesh_color.row(0) = Eigen::RowVector3d(255, 255, 120) / 255.0;
//	reset_viewer();
//
//	// Show points and align camera
//	points_data().clear();
//	points_data().set_points(state.mesh.points, Eigen::RowVector3d(1, 0, 0));
//	viewer.core.align_camera_center(state.mesh.points);
//	double extent = (state.mesh.points.colwise().maxCoeff() - state.mesh.points.colwise().minCoeff()).maxCoeff();
//	points_data().point_size = float(0.008 * extent);
//	fix_color(points_data());
//
//	// Compute and show convex hull + triangulation
//	compute_hull();
//	//clean_hull();
//	compute_triangulation();
//	viewer_control();
//	return true;
//}

bool UIState::load() {
	//assert(false);
	if (!state.load(save_dir)) { return false; }
	current_region_status = "";
	//img.resize(0, 0);
	// reset flags

	mesh_color.resize(1, 3);
	mesh_color.row(0) = Eigen::RowVector3d(255, 255, 120) / 255.0;
	reset_viewer();

	color_code = true;

	selected_region = -1;
	viewer_control();
	return true;
}

void UIState::detect_vertices() {
	state.detect_vertices();

	mesh_color.resize(1, 3);
	mesh_color.row(0) = Eigen::RowVector3d(255, 255, 120) / 255.0;
	//reset_viewer();
	show_points = true;

	compute_hull();
	//clean_hull();
	compute_triangulation();

	color_code = true;
	selected_region = -1;
	current_region_status = "";
	viewer_control();
}

//bool UIState::load_param(std::string name)
//{
//	return state.load_param(name);
//}

bool UIState::save() {
	int nError = 0;
	cellogram_mkdir(save_dir);
	cellogram_mkdir(save_dir + "/cellogram");
	return state.save(save_dir + "/cellogram");
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
	vertex_color = Eigen::RowVector3f(1,0,0);
	show_hull = true;
	show_points = true;
	show_mesh_fill = true;
	show_image = true;
	show_matching = false;
	show_bad_regions = false;
	show_mesh_fill = false;
}

//TODO refactor when more clear
void UIState::load_image(std::string fname) {
	
	state.load_image(fname);

	const int index = fname.find_last_of(".");
	save_dir = fname.substr(0, index);

	texture = (state.img.array() * 255).cast<unsigned char>();

	selected_region = -1;
	current_region_status = "";

	color_code = false;
	show_mesh_fill = false;
	show_image = true;

	int xMax = state.img.cols();
	int yMax = state.img.rows();
	Eigen::MatrixXd V(4, 3);
	V <<
		0, 0, 0,
		yMax, 0, 0,
		yMax, xMax, 0,
		0, xMax, 0;

	viewer.core.align_camera_center(V);

	double extent = (V.colwise().maxCoeff() -V.colwise().minCoeff()).maxCoeff();
	//HIGH dpi

	//int width, height;
	//glfwGetFramebufferSize(viewer.window, &width, &height);
	//int width_window, height_window;
	//glfwGetWindowSize(viewer.window, &width_window, &height_window);
	//const int highdpi = width / width_window;


	std::cout << extent << std::endl;
	points_data().point_size =  float(700. / extent)+5;


	viewer_control();
}

void UIState::display_image()
{
	int xMax = state.img.cols();
	int yMax = state.img.rows();

	// Replace the mesh with a triangulated square
	img_V.resize(4, 3);
	img_V <<
		0, 0, 0,
		yMax, 0, 0,
		yMax, xMax, 0,
		0, xMax, 0;
	img_F.resize(2, 3);
	img_F <<
		0, 1, 2,
		2, 3, 0;
	Eigen::MatrixXd UV(4, 2);
	UV <<
		0, 0,
		1, 0,
		1, 1,
		0, 1;

	image_data().set_mesh(img_V, img_F);
	image_data().set_uv(UV);
	image_data().show_faces = true;
	image_data().show_lines = false;
	image_data().show_texture = true;
	// Use the image as a texture
	image_data().set_texture(texture, texture, texture);
	image_data().set_colors(Eigen::RowVector3d(1, 1, 1));
}

void UIState::viewer_control()
{
	hull_data().clear();
	points_data().clear();
	image_data().clear();
	matching_data().clear();
	bad_region_data().clear();
	selected_data().clear();

	// Read all the UIState flags and update display
	// hull

	int n = (int)state.hull_polygon.rows();

	if (show_hull && n >= 3) {
		// Draw edges

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
	Eigen::MatrixXd V = t * state.mesh.moved + (1 - t) * state.mesh.moved;

	if (V.size() > 0)
		points_data().set_mesh(V, state.mesh.triangles);

	points_data().show_lines = dragging_id < 0;

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
		//if (state.regions.empty())
		//{
		//	Eigen::MatrixXd C(V.rows(), 3);
		//	C.col(0).setConstant(vertex_color(0));
		//	C.col(1).setConstant(vertex_color(1));
		//	C.col(2).setConstant(vertex_color(2));

		//	if (color_code)
		//	{
		//		Eigen::VectorXd param(V.rows());
		//		for (int i = 0; i < V.rows(); i++)
		//		{
		//			param(i) = state.mesh.params.std_x(i);
		//		}
		//		//igl::parula(param, true, C);
		//		igl::ColorMapType cm = igl::ColorMapType::COLOR_MAP_TYPE_INFERNO;
		//		igl::colormap(cm, param, true, C);
		//	}

		//	points_data().set_points(V, C);
		//}
		//else
		//{
			Eigen::MatrixXd C(V.rows(), 3);
			C.col(0).setConstant(vertex_color(0));
			C.col(1).setConstant(vertex_color(1));
			C.col(2).setConstant(vertex_color(2));

			if (color_code)
			{
				Eigen::VectorXd param(V.rows());
				for (int i = 0; i < V.rows(); i++)
				{
					param(i) = state.mesh.params.pval_Ar(i);
				}

				igl::ColorMapType cm = igl::ColorMapType::COLOR_MAP_TYPE_INFERNO;
				igl::colormap(cm,param, true, C);
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

		//}
	}

	// fill
	points_data().show_faces = show_mesh_fill;
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

	// show selected region
	if (show_selected_region && selected_region > 0)
	{
		selected_data().clear();
		int nTri = state.regions[selected_region].region_triangles.size();
		Eigen::MatrixXi region_tri(nTri, 3);
		for (int j = 0; j < nTri; ++j)
		{
			region_tri.row(j) = state.mesh.triangles.row(state.regions[selected_region].region_triangles(j));
		}
		//Eigen::MatrixXd vtmp = V;
		//vtmp.col(2).setConstant(0.001);
		Eigen::RowVector3d color;
		color << 231, 76, 60;
		color /= 255;
		//selected_data().set_mesh(V, region_tri);
		//selected_data().set_colors(color);
		//selected_data().show_faces = true;

		int nEdgePts = state.regions[selected_region].region_boundary.size();
		Eigen::MatrixXd region_edge1(nEdgePts,3), region_edge2(nEdgePts, 3);
		for (int j = 0; j < nEdgePts; ++j)
		{
			region_edge1.row(j) = V.row(state.regions[selected_region].region_boundary(j));
			region_edge2.row(j) = V.row(state.regions[selected_region].region_boundary((j + 1) % nEdgePts));
		}

		selected_data().add_edges(region_edge1, region_edge2, color);
		selected_data().line_width = 3.0;
	}

	// Fix shininess for all layers
	fix_color(hull_data());
	fix_color(points_data());
	fix_color(image_data());
	fix_color(bad_region_data());
	fix_color(matching_data());
	fix_color(selected_data());
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


void UIState::create_region_label()
{
	//mesh_color.resize(state.mesh.points.rows(),3);
	mesh_color.resize(state.mesh.triangles.rows(), 3);
	mesh_color.setOnes();

	for (int i = 0; i < state.regions.size(); ++i)
	{
		Eigen::RowVector3d color;

		switch (state.regions[i].status)
		{
		case Region::OK:
			color << 46, 204, 113; break;
		case Region::TOO_MANY_VERTICES:
			color << 155, 89, 182; break;
		case Region::TOO_FEW_VERTICES:
			color << 241, 196, 15; break;
		case Region::REGION_TOO_LARGE:
			color << 41, 128, 185; break;
		case Region::NO_SOLUTION:
			color << 192, 57, 43; break;
		case Region::NOT_PROPERLY_CLOSED:
			color << 149, 165, 166; break;
		default:
			color << 52, 73, 94;  break;
		}
		color /= 255;
		//for (int j = 0; j < state.regions[i].region_interior.size(); ++j)
		//{
		//	mesh_color.row(state.regions[i].region_interior(j)) = color;
		//}
		for (int j = 0; j < state.regions[i].region_triangles.size(); ++j)
		{
			mesh_color.row(state.regions[i].region_triangles(j)) = color;
		}
	}

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
