#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "cellogram/State.h"

#include <Eigen/Dense>

#include <igl/colormap.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

// -----------------------------------------------------------------------------

class UIState : public igl::opengl::glfw::imgui::ImGuiMenu {
private:
	UIState();

public:
	typedef igl::opengl::glfw::imgui::ImGuiMenu super;
	static UIState &ui_state();

	virtual ~UIState() = default;

public:
	// Viewer
	igl::opengl::glfw::Viewer viewer;

	// Scene
	State &state;

	// Multiple meshes ids
	int hull_id;
	int points_id;
	int image_id;
	int bad_region_id;
	int matching_id;
	int selected_id;
	int physical_id;


	int selected_region = -1;
	int selected_param = 0;
	int dragging_id = -1;
	// UI options
	// double foo;
	
	Eigen::MatrixXi img_F;
	Eigen::MatrixXd img_V;
	Eigen::MatrixXf hist;

	// Display flags
	float t;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture;
	Eigen::RowVector3f vertex_color;
	Eigen::MatrixXd mesh_color;

	bool image_from_pillars = false;

	bool show_mesh = true;
	bool show_hull = false;
	//bool image_loaded = false;
	bool show_points = false;
	bool show_mesh_fill = true;
	bool show_image = false;
	bool show_matching = false;
	bool show_bad_regions = false;
	bool color_code = false;
	bool show_selected_region = true;
	bool analysis_mode = false;
	bool show_traction_forces = true;
	
	// 3d visualizer
	int view_mode_3d = 0;
	static const int NO_VIEW_SELECTED = 0;
	static const int X_DISP_SELECTED = 1;
	static const int Y_DISP_SELECTED = 2;
	static const int Z_DISP_SELECTED = 3;
	static const int MAG_DISP_SELECTED = 4;

	// Clicking flags
	bool select_region = false;
	bool add_vertex = false;
	bool delete_vertex = false;
	bool make_vertex_good = false;
	bool make_vertex_bad = false;
	bool move_vertex = false;
	int split_region = -1;
	Eigen::Vector2i split_end_points;
	// Image
	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> img;

	std::string current_region_status;

	std::string save_dir = "";
	bool data_available = false;
public:
	void initialize();

	void launch();

	igl::opengl::ViewerData & mesh_by_id(int id);

	bool load();
	//bool load_param(std::string name);

	bool save();

	virtual bool mouse_down(int button, int modifier) override;
	virtual bool mouse_move(int button, int modifier) override;
	virtual bool mouse_up(int button, int modifier) override;
	virtual bool mouse_scroll(float delta_y) override;

	void load_image(std::string name);
	void detect_vertices();
	void display_image();
	void compute_hull();
	void clean_hull();
	void compute_triangulation();
	void compute_histogram();

	void export_region();
	void reset_viewer();
	void deselect_all_buttons();
public:
	igl::opengl::ViewerData & points_data() { return mesh_by_id(points_id); }
	igl::opengl::ViewerData & hull_data() { return mesh_by_id(hull_id); }
	igl::opengl::ViewerData & image_data() { return mesh_by_id(image_id); }
	igl::opengl::ViewerData & bad_region_data() { return mesh_by_id(bad_region_id); }
	igl::opengl::ViewerData & matching_data() { return mesh_by_id(matching_id); }
	igl::opengl::ViewerData & selected_data() { return mesh_by_id(selected_id); }
	igl::opengl::ViewerData & physical_data() { return mesh_by_id(physical_id); }

private:
	bool block_mouse_behavior(int button);
	void viewer_control();
	void viewer_control_2d();
	void viewer_control_3d();
	void draw_mesh();
	void fix_color(igl::opengl::ViewerData &data);
	void create_region_label();
	void build_region_edges(const Eigen::MatrixXd &pts, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2, Eigen::MatrixXd &C);

	void build_menu_bar();

	bool show_file_menu = true;
	void draw_file_menu(int x, int y, int &y_return);

	bool show_points_menu = true;
	void draw_points_menu(int x, int y);

	bool show_mesh_menu = true;
	void draw_mesh_menu(int x, int y);

	bool show_analysis_menu = true;
	void draw_analysis_menu(int x, int y);

	bool show_histogram = true;
	void draw_histogram(int x, int y);

	bool show_legend = true;
	void draw_legend(int x, int y);

	bool show_view_options = true;
	void draw_view_options(int x, int y);

	bool show_region_options = true;
	void draw_region_menu(int x, int y);
public:
	// Menu stuff
	void draw_viewer_window() override { }
	void draw_viewer_menu() override;
	void draw_custom_window() override;

	bool key_pressed(unsigned int unicode_key, int modifiers) override;
	bool key_up(int key, int modifiers) override;
};

// -----------------------------------------------------------------------------

} // namespace cellogram
