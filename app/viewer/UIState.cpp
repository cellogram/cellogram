////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include "FileDialog.h"
#include <cellogram/PolygonUtils.h>
#include <cellogram/StringUtils.h>
#include <igl/colon.h>
#include <igl/colormap.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <imgui/imgui.h>
#include <GLFW/glfw3.h>
#include <sys/stat.h>  // no clue why required -- man pages say so
#include <sys/types.h> // required for stat.h
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace {

		int cellogram_mkdir(const std::string &path) {
			int nError;
	#if defined(_WIN32)
			std::wstring widestr = std::wstring(path.begin(), path.end());
		nError = _wmkdir(widestr.c_str()); // can be used on Windows
	#else
		nError = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // can be used on non-Windows
	#endif
		return nError;
	}

	void get_region_color(const int &status, Eigen::RowVector3d &color) {
		switch (status) {
			case Region::OK:
			color << 46, 204, 113;
			break;
			case Region::TOO_MANY_VERTICES:
			color << 155, 89, 182;
			break;
			case Region::TOO_FEW_VERTICES:
			color << 241, 196, 15;
			break;
			case Region::REGION_TOO_LARGE:
			color << 41, 128, 185;
			break;
			case Region::NO_SOLUTION:
			color << 192, 57, 43;
			break;
			case Region::NOT_PROPERLY_CLOSED:
			color << 149, 165, 166;
			break;
			default:
			color << 52, 73, 94;
			break;
		}
		color /= 255;
	}

} // anonymous namespace

// -----------------------------------------------------------------------------

UIState::UIState() : state(State::state()) { reset_viewer(); }

UIState &UIState::ui_state() {
	static UIState instance;
	return instance;
}

bool UIState::mouse_move(int button, int modifier) {
	if (dragging_id == -1)
		return false;

	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;

	// reset update position with x and y, store prevx and y to compute delatax,
	// (use current zoom to move point accordingly???)
	int fid;
	Eigen::Vector3f bc;
	igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model, viewer.core.proj,
		viewer.core.viewport, img_V, img_F, fid, bc);

	double xNew = 0, yNew = 0, zNew = 0;
	for (int i = 0; i < 3; i++) {
		// xNew += bc(i) * V(state.mesh.triangles(fid, i), 0);
		xNew += bc(i) * img_V(img_F(fid, i), 0);
		yNew += bc(i) * img_V(img_F(fid, i), 1);
		zNew += 0;
	}

	state.mesh.moved.row(dragging_id) = Eigen::RowVector3d(xNew, yNew, zNew);

	viewer_control();

	return true;
}

bool UIState::block_mouse_behavior(int button) {
	if (analysis_mode)
		return false;

	if (button == 0) {
		this->viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::None;
		return true;
	}

	return false;
}

bool UIState::mouse_up(int button, int modifier) {
	if (button != 0 || dragging_id == -1)
		return false;
	// move_vertex = false;
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

	// if (!select_region && !add_vertex && !delete_vertex && !make_vertex_good &&
	// !make_vertex_bad) {

	//	if (button == 0)
	//	{
	//		this->viewer.mouse_mode =
	// igl::opengl::glfw::Viewer::MouseMode::None; 		return true;
	//	}

	//	return false;
	//}

	int fid;
	Eigen::Vector3f bc;
	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;
	Eigen::MatrixXd V = t * state.mesh.points + (1 - t) * state.mesh.moved;

	if (V.size() <= 0)
		return block_mouse_behavior(button);

	if (select_region) {
		if (!igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model, viewer.core.proj,
			viewer.core.viewport, V, state.mesh.triangles, fid, bc))
		{
			return block_mouse_behavior(button);
		}

		select_region = false;

		// set selected_region
		for (int i = 0; i < state.regions.size(); i++) {
			for (int j = 0; j < state.regions[i].region_triangles.rows(); j++) {
				if (fid == state.regions[i].region_triangles(j)) {
					selected_region = i;
					viewer_control();
					current_region_status = Region::pretty_status(state.regions[i].status);
					return true;
				}
			}
		}
	} else {
		if (!igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model, viewer.core.proj,
			viewer.core.viewport, img_V, img_F, fid, bc))
		{
			return block_mouse_behavior(button);
		}

		double xNew = 0, yNew = 0, zNew = 0;
		for (int i = 0; i < 3; i++) {
			// xNew += bc(i) * V(state.mesh.triangles(fid, i), 0);
			xNew += bc(i) * img_V(img_F(fid, i), 0);
			yNew += bc(i) * img_V(img_F(fid, i), 1);
			zNew += 0;
		}

		Eigen::VectorXd delta;
		delta = (V.col(0).array() - xNew).cwiseAbs2() + (V.col(1).array() - yNew).cwiseAbs2();
		int vid;
		delta.minCoeff(&vid);

		// add_vertex = false;
		if (add_vertex) {
			state.add_vertex(Eigen::Vector3d(xNew, yNew, zNew));

			state.reset_state();
			reset_viewer();
			viewer_control();
			return true;
		} else if (delete_vertex) {
			// find vertex in V closest to xNew,yNew
			state.delete_vertex(vid);

			state.reset_state();
			reset_viewer();
			viewer_control();
			return true;
		} else if (move_vertex) {
			dragging_id = vid;
			return true;
		} else if (make_vertex_good || make_vertex_bad) {
			// find maximum barycenter coordinates and get vertex id
			if (make_vertex_good) {
				state.mesh.vertex_status_fixed(vid) = 1;
			} else if (make_vertex_bad) {
				state.mesh.vertex_status_fixed(vid) = -1;
			}

			make_vertex_good = false;
			make_vertex_bad = false;

			state.detect_bad_regions();
			state.check_regions();
			create_region_label();
			viewer_control();

			return true;
		} else if (split_region > -1) {
			split_end_points(split_region) = vid;
			split_region++;

			viewer_control();

			if (split_region > 1) {
				state.split_region(split_end_points);
				split_region = -1;
				viewer_control();
			}
			return true;
		}
	}

	return block_mouse_behavior(button);
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
	viewer.resize(1400, 1000);
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
	viewer.append_mesh();
	viewer.append_mesh();
	viewer.data_list[0].id = hull_id = 0;
	viewer.data_list[1].id = points_id = 1;
	viewer.data_list[2].id = image_id = 2;
	viewer.data_list[3].id = bad_region_id = 3;
	viewer.data_list[4].id = matching_id = 4;
	viewer.data_list[5].id = selected_id = 5;
	viewer.data_list[6].id = physical_id = 6;

	assert(viewer.data_list.size() == physical_id + 1);
}

void UIState::init(igl::opengl::glfw::Viewer *_viewer) {
	super::init(_viewer);

	ImGuiIO& io = ImGui::GetIO();
	io.IniFilename = nullptr;

	glfwSetWindowTitle(viewer.window, "Cellogram viewer");
}

void UIState::launch() {
	// Launch viewer
	viewer.launch();
}

////////////////////////////////////////////////////////////////////////////////

igl::opengl::ViewerData &UIState::mesh_by_id(int id) {
	size_t index = viewer.mesh_index(id);
	assert(viewer.data_list[index].id == id);
	return viewer.data_list[index];
}

// bool UIState::load(std::string name) {
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
//	points_data().set_points(state.mesh.points, Eigen::RowVector3d(1, 0,
// 0)); 	viewer.core.align_camera_center(state.mesh.points); 	double
// extent = (state.mesh.points.colwise().maxCoeff() -
// state.mesh.points.colwise().minCoeff()).maxCoeff(); points_data().point_size
// = float(0.008 * extent); 	fix_color(points_data());
//
//	// Compute and show convex hull + triangulation
//	compute_hull();
//	//clean_hull();
//	compute_triangulation();
//	viewer_control();
//	return true;
//}

bool UIState::load() {
	// assert(false);
	if (!state.load(save_dir)) {
		return false;
	}
	current_region_status = "";
	// img.resize(0, 0);
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

	if (state.mesh.points.size() == 0)
		return;

	mesh_color.resize(1, 3);
	mesh_color.row(0) = Eigen::RowVector3d(255, 255, 120) / 255.0;
	// reset_viewer();
	show_points = true;

	compute_hull();
	// clean_hull();
	compute_triangulation();

	color_code = true;
	selected_region = -1;
	current_region_status = "";
	viewer_control();
}

// bool UIState::load_param(std::string name)
//{
//	return state.load_param(name);
//}

bool UIState::save() {
	int nError = 0;
	cellogram_mkdir(save_dir);
	cellogram_mkdir(save_dir + "/cellogram");
	bool ok = state.save(save_dir + "/cellogram");
	data_available = state.is_data_available(save_dir);

	return ok;
}

bool UIState::mouse_scroll(float delta_y) {
	viewer.scroll_position += delta_y;

	// Only zoom if there's actually a change
	if (delta_y != 0) {
		float mult = (1.0 + ((delta_y > 0) ? 1. : -1.) * 0.1);
		const float min_zoom = 0.1f;
		viewer.core.camera_zoom =
		(viewer.core.camera_zoom * mult > min_zoom ? viewer.core.camera_zoom * mult : min_zoom);
	}

	return super::mouse_scroll(delta_y);
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

bool UIState::key_pressed(unsigned int unicode_key, int modifiers) {
	if (modifiers != 0) {
		return super::key_pressed(unicode_key, modifiers);
	}

	switch (unicode_key) {
		case 'a':
		case 'A':
		delete_vertex = false;
		add_vertex = !add_vertex;
		return true;
		case 'd':
		case 'D':
		add_vertex = false;
		delete_vertex = !delete_vertex;
		return true;
	}

	return super::key_pressed(unicode_key, modifiers);
}

bool UIState::key_up(int key, int modifiers) {
	if (modifiers == 2) {
		if (key == 'o' || key == 'O') {
			std::string fname = FileDialog::openFileName(DATA_DIR, {"*.png", "*.tif", "*.tiff"});
			if (!fname.empty()) {
				load_image(fname);

				show_image = true;

				// update UI
				viewer_control();
			}

			return true;
		}
		if (key == 's' || key == 'S') {
			save();

			return true;
		}
	}

	return super::key_up(key, modifiers);
}
// -----------------------------------------------------------------------------

void UIState::compute_triangulation() {
	// State
	// state.compute_triangulation();

	// UI
	show_points = true;
}

void UIState::compute_histogram() {
	int n_bins = 64;
	hist.resize(n_bins + 1, 1);
	hist.setZero();

	for (long i = 0; i < state.img.size(); ++i) {
		hist(state.img(i) * n_bins)++;
	}
}

void UIState::export_region() {
	// Eigen::VectorXi boundary = state.regions[selected_region].region_boundary;
	// Eigen::VectorXi internal = state.regions[selected_region].region_interior;
	//
	// Eigen::MatrixXd V(boundary.size() + internal.size(), 3);
	// for (int i = 0; i < boundary.size(); i++)
	//{
	//	V(i, 0) = state.mesh.detected(boundary(i), 0);
	//	V(i, 1) = state.mesh.detected(boundary(i), 1);
	//	double std_x = state.mesh.params.std_x(boundary(i));
	//	double std_y = state.mesh.params.std_y(boundary(i));
	//	V(i, 2) = std::sqrt(std_x*std_x + std_y * std_y);
	//}
	// for (int i = 0; i < internal.size(); i++)
	//{
	//	V(i + boundary.size(), 0) = state.mesh.detected(internal(i), 0);
	//	V(i + boundary.size(), 1) = state.mesh.detected(internal(i), 1);
	//	double std_x = state.mesh.params.std_x(internal(i));
	//	double std_y = state.mesh.params.std_y(internal(i));
	//	V(i + boundary.size(), 2) = std::sqrt(std_x*std_x + std_y*std_y);
	//}

	//// find vertices whose edges can be regarded as possibly wrong
	// Eigen::VectorXi wrong_boundary =
	// state.increase_boundary(state.mesh.boundary); wrong_boundary =
	// state.increase_boundary(wrong_boundary);

	// int index = 0;
	// Eigen::MatrixXi E(boundary.size(), 2);
	// for (int i = 0; i < boundary.size(); i++)
	//{
	//	if ((wrong_boundary.array() - boundary(i)).abs().minCoeff() == 0 ||
	//(wrong_boundary.array() - boundary((i + 1) %
	// boundary.size())).abs().minCoeff() == 0) 		continue; 	E(index, 0) =
	// i; 	E(index, 1) = (i + 1) % boundary.size(); 	index++;
	//}
	// E.conservativeResize(index, 2);

	// cellogram_mkdir(save_dir);
	// cellogram_mkdir(save_dir + "/cellogram");
	// cellogram_mkdir(save_dir + "/cellogram" + "/regions");
	// cellogram_mkdir(save_dir + "/cellogram" + "/regions/" +
	// std::to_string(selected_region));

	// std::string path = save_dir + "/cellogram" + "/regions/" +
	// std::to_string(selected_region);
	//
	//{
	//	std::ofstream V_path(path + "/V.vert");
	//	V_path << V << std::endl;
	//	V_path.close();
	//}
	//{
	//	std::ofstream E_path(path + "/E.edge");
	//	E_path << E << std::endl;
	//	E_path.close();
	//}

	////---------- Save entire image and vector<bool> for faces that may not be
	/// changed
	cellogram_mkdir(save_dir);
	cellogram_mkdir(save_dir + "/cellogram");
	cellogram_mkdir(save_dir + "/cellogram" + "/untangler");
	std::string path = save_dir + "/cellogram/untangler";

	{
		std::ofstream V_path(path + "/V.vert");
		V_path << state.mesh.moved << std::endl;
		V_path.close();
	}
	{
		std::ofstream F_path(path + "/F.tri");
		F_path << state.mesh.triangles << std::endl;
		F_path.close();
	}

	Eigen::VectorXi fixed_face(state.mesh.triangles.rows());
	fixed_face.setOnes();

	// Loop through regions and set them to 0
	for (auto &r : state.regions) {
		for (int i = 0; i < r.region_triangles.size(); i++) {
			fixed_face(r.region_triangles(i)) = 0;
		}
	}

	// Loop through boundary and also set them to zero
	Eigen::VectorXi boundary = state.increase_boundary(state.mesh.boundary);
	// boundary = state.increase_boundary(boundary);
	for (int i = 0; i < boundary.size(); i++) {
		std::vector<int> fixed_ind = state.mesh.vertex_to_tri[boundary(i)];
		for (int j = 0; j < (int) fixed_ind.size(); j++) {
			fixed_face(fixed_ind[j]) = 0;
		}
	}
	{
		std::ofstream fixed_path(path + "/fixed.txt");
		fixed_path << fixed_face << std::endl;
		fixed_path.close();
	}
}

void UIState::reset_viewer() {
	// Display flags
	t = 0;
	vertex_color = Eigen::RowVector3f(1, 0, 0);
	selected_region = -1;
	show_hull = true;
	show_points = true;
	show_mesh_fill = true;
	show_image = true;
	show_matching = false;
	show_bad_regions = false;
	show_mesh_fill = false;
}

void UIState::deselect_all_buttons() {
	select_region = false;
	add_vertex = false;
	delete_vertex = false;
	make_vertex_good = false;
	make_vertex_bad = false;
	move_vertex = false;
}

// TODO refactor when more clear
void UIState::load_image(std::string fname) {
	state.load_image(fname);

	const int index = fname.find_last_of(".");
	save_dir = fname.substr(0, index);
	data_available = state.is_data_available(save_dir);

	texture = (state.img.array() * 255).cast<unsigned char>();

	selected_region = -1;
	current_region_status = "";

	color_code = false;
	analysis_mode = false;
	show_mesh_fill = false;
	show_image = true;

	int xMax = state.img.cols();
	int yMax = state.img.rows();
	Eigen::MatrixXd V(4, 3);
	V << 0, 0, 0, yMax, 0, 0, yMax, xMax, 0, 0, xMax, 0;

	viewer.core.align_camera_center(V);
	compute_histogram();

	double extent = (V.colwise().maxCoeff() - V.colwise().minCoeff()).maxCoeff();
	points_data().point_size = std::ceil(float(700. / extent)) + 3;

	// HIGH dpi
	// int width, height;
	// glfwGetFramebufferSize(viewer.window, &width, &height);
	// int width_window, height_window;
	// glfwGetWindowSize(viewer.window, &width_window, &height_window);
	// const int highdpi = width / width_window;

	viewer_control();
}

void UIState::display_image() {
	int xMax = state.img.cols();
	int yMax = state.img.rows();
	const float offset = 0.5;
	// Replace the mesh with a triangulated square
	img_V.resize(4, 3);
	img_V << offset, offset, 0, yMax + offset, offset, 0, yMax + offset, xMax + offset, 0, offset, xMax + offset, 0;

	img_F.resize(2, 3);
	img_F << 0, 1, 2, 2, 3, 0;

	Eigen::MatrixXd UV(4, 2);
	UV << 0, 1, 1, 1, 1, 0, 0, 0;

	image_data().set_mesh(img_V, img_F);
	image_data().set_uv(UV);
	image_data().show_faces = true;
	image_data().show_lines = false;
	image_data().show_texture = true;
	// Use the image as a texture
	image_data().set_texture(texture, texture, texture);
	image_data().set_colors(Eigen::RowVector3d(1, 1, 1));
}

void UIState::viewer_control() {
	if (analysis_mode) {
		viewer_control_3d();
		viewer.core.orthographic = false;
	} else {
		viewer_control_2d();
		viewer.core.trackball_angle = Eigen::Quaternionf::Identity();
		viewer.core.orthographic = true;
	}
}

void UIState::viewer_control_2d() {
	viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_NO_ROTATION);

	hull_data().clear();
	points_data().clear();
	image_data().clear();
	matching_data().clear();
	bad_region_data().clear();
	selected_data().clear();
	physical_data().clear();

	// Read all the UIState flags and update display
	// hull

	int n = (int)state.hull_polygon.rows();

	if (show_hull && n >= 3) {
		// Draw edges

		Eigen::MatrixXd P2;
		P2.resizeLike(state.hull_polygon);
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
	Eigen::MatrixXd V = t * state.mesh.points + (1 - t) * state.mesh.moved;

	if (V.size() > 0 && show_mesh)
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
	if (show_points) {
		if (viewer.window != NULL) {
			float ui_scaling_factor = hidpi_scaling() / pixel_ratio();
			points_data().point_size = std::ceil(ui_scaling_factor) + 3;
		}

		// if (state.regions.empty())
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
		//		igl::ColorMapType cm =
		// igl::ColorMapType::COLOR_MAP_TYPE_INFERNO; 		igl::colormap(cm, param,
		// true, C);
		//	}

		//	points_data().set_points(V, C);
		//}
		// else
		//{
		Eigen::MatrixXd C(V.rows(), 3);
		C.col(0).setConstant(vertex_color(0));
		C.col(1).setConstant(vertex_color(1));
		C.col(2).setConstant(vertex_color(2));
		if (color_code) {
			Eigen::VectorXd param(V.rows());
			for (int i = 0; i < V.rows(); i++) {
				// param(i) = state.mesh.params.pval_Ar(i);
				param(i) = std::sqrt(state.mesh.params.std_x(i) * state.mesh.params.std_x(i) +
					state.mesh.params.std_y(i) * state.mesh.params.std_y(i));
			}

			igl::ColorMapType cm = igl::ColorMapType::COLOR_MAP_TYPE_INFERNO;
			igl::colormap(cm, param, true, C);
		} else {
			for (auto &r : state.regions) {
				for (int i = 0; i < r.region_boundary.size(); ++i) {
					C.row(r.region_boundary(i)) = Eigen::RowVector3d(0, 0, 1);
				}

				for (int i = 0; i < r.region_interior.size(); ++i) {
					C.row(r.region_interior(i)) = Eigen::RowVector3d(0, 1, 0);
				}
			}
		}

		points_data().set_points(V, C);

		//}
	}

	// fill
	points_data().show_faces = show_mesh_fill;
	// if (show_mesh_fill && !state.regions.empty())
	if (show_mesh_fill)
		create_region_label();
	points_data().set_colors(mesh_color);

	// bad regions
	if (show_bad_regions) {
		Eigen::MatrixXd bad_P1, bad_P2, C;
		build_region_edges(V, bad_P1, bad_P2, C);
		if (show_mesh_fill) {
			bad_region_data().add_edges(bad_P1, bad_P2, Eigen::RowVector3d(0, 0, 0));
			bad_region_data().line_width = 3.0;
		} else {
			bad_region_data().add_edges(bad_P1, bad_P2, C);
			bad_region_data().line_width = 2.0;
		}
	}

	// image
	if (show_image) {
		if (state.img.size() <= 0) {
			show_image = false;
		} else {
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
	if (show_selected_region && selected_region != -1) {
		// selected_data().clear();
		// int nTri = state.regions[selected_region].region_triangles.size();
		// Eigen::MatrixXi region_tri(nTri, 3);
		// for (int j = 0; j < nTri; ++j)
		//{
		//	region_tri.row(j) =
		// state.mesh.triangles.row(state.regions[selected_region].region_triangles(j));
		//}
		// Eigen::MatrixXd vtmp = V;
		// vtmp.col(2).setConstant(0.001);
		Eigen::RowVector3d color;
		color << 231, 76, 60;
		color /= 255;
		// selected_data().set_mesh(V, region_tri);
		// selected_data().set_colors(color);
		// selected_data().show_faces = true;

		int nEdgePts = state.regions[selected_region].region_boundary.size();
		Eigen::MatrixXd region_edge1(nEdgePts, 3), region_edge2(nEdgePts, 3);
		for (int j = 0; j < nEdgePts; ++j) {
			region_edge1.row(j) = V.row(state.regions[selected_region].region_boundary(j));
			region_edge2.row(j) = V.row(state.regions[selected_region].region_boundary((j + 1) % nEdgePts));
		}

		region_edge1.col(2).array() += 1e-2;
		region_edge2.col(2).array() += 1e-2;

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

void UIState::viewer_control_3d() {
	viewer_control_2d();
	physical_data().clear();

	// ROTATION_TYPE_TRACKBALL
	// ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP
	viewer.core.set_rotation_type(igl::opengl::ViewerCore::RotationType::ROTATION_TYPE_TRACKBALL);

	if (state.mesh3d.V.size() == 0)
		return;

	Eigen::MatrixXd Vtmp = state.mesh3d.V / state.mesh.scaling;

	// Vtmp.col(2).array() -= 0.1;
	// physical_data().set_mesh(Vtmp, state.mesh3d.F);


	Eigen::MatrixXd C, disp;

	const auto &fun = show_traction_forces ? state.mesh3d.traction_forces : state.mesh3d.displacement;

	switch (view_mode_3d) {
		case Mesh3DAttribute::NONE:
		{
			C = Eigen::RowVector3d(129. / 255, 236. / 255, 236. / 255);
			min_val = max_val = 0;
			break;
		}
	 	case Mesh3DAttribute::X_DISP:
		{
			const auto v_fun = fun.col(0).eval();
			min_val = v_fun.minCoeff();
			max_val = v_fun.maxCoeff();
			igl::colormap(cm, v_fun, true, C);
			break;
		}
		case Mesh3DAttribute::Y_DISP:
		{
			const auto v_fun = fun.col(1).eval();
			min_val = v_fun.minCoeff();
			max_val = v_fun.maxCoeff();
			igl::colormap(cm, v_fun, true, C);
			break;
		}
		case Mesh3DAttribute::Z_DISP:
		{
			const auto v_fun = fun.col(2).eval();
			min_val = v_fun.minCoeff();
			max_val = v_fun.maxCoeff();
			igl::colormap(cm, v_fun, true, C);
			break;
		}
		case Mesh3DAttribute::NORM_DISP:
		{
			const auto v_fun = fun.rowwise().norm().eval();
			min_val = v_fun.minCoeff();
			max_val = v_fun.maxCoeff();
			igl::colormap(cm, v_fun, true, C);
			break;
		}
	}

	// std::cout << C << std::endl;
	// std::cout << "\n\ncol\n" << fun.col(0).eval() << std::endl;

	Eigen::MatrixXd normals;

	if (state.image_from_pillars) {
		physical_data().set_points(Vtmp, C);
	} else {
		Vtmp.col(2).array() -= 0.1;
		physical_data().set_mesh(Vtmp, state.mesh3d.F);

		igl::per_face_normals(Vtmp, state.mesh3d.F, normals);
		// normals *= -1;
		physical_data().set_normals(normals);
		physical_data().set_colors(C);
	}

	physical_data().show_lines = view_mode_3d == Mesh3DAttribute::NONE;
	// fix_color(physical_data());
}

void UIState::draw_mesh() {
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

void UIState::create_region_label() {
	// mesh_color.resize(state.mesh.points.rows(),3);
	mesh_color.resize(state.mesh.triangles.rows(), 3);
	mesh_color.setOnes();

	for (int i = 0; i < (int) state.regions.size(); ++i) {
		Eigen::RowVector3d color;
		get_region_color(state.regions[i].status, color);

		// for (int j = 0; j < state.regions[i].region_interior.size(); ++j)
		//{
		//	mesh_color.row(state.regions[i].region_interior(j)) = color;
		//}
		for (int j = 0; j < state.regions[i].region_triangles.size(); ++j) {
			mesh_color.row(state.regions[i].region_triangles(j)) = color;
		}
	}
}

void UIState::build_region_edges(const Eigen::MatrixXd &pts, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2,
	Eigen::MatrixXd &C)
{
	bad_P1.resize(pts.rows(), 3);
	bad_P2.resize(pts.rows(), 3);
	C.setZero(pts.rows(), 3);

	int index = 0;
	int c_index = 0;
	Eigen::MatrixXd local_bad_P1, local_bad_P2;
	for (auto &r : state.regions) {
		r.compute_edges(pts, local_bad_P1, local_bad_P2);
		bad_P1.block(index, 0, local_bad_P1.rows(), local_bad_P1.cols()) = local_bad_P1;
		bad_P2.block(index, 0, local_bad_P2.rows(), local_bad_P2.cols()) = local_bad_P2;

		Eigen::RowVector3d color;
		get_region_color(r.status, color);

		for (int j = 0; j < r.region_boundary.size(); j++) {
			C.row(c_index) = color;
			c_index++;
		}

		index += local_bad_P2.rows();
	}
	bad_P1.conservativeResize(index, 3);
	bad_P2.conservativeResize(index, 3);
}

// -----------------------------------------------------------------------------

} // namespace cellogram
