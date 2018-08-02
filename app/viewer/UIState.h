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
	enum class Mesh3DAttribute : int {
		NONE,
		X_DISP,
		Y_DISP,
		Z_DISP,
		NORM_DISP,
	};
	Mesh3DAttribute view_mode_3d = Mesh3DAttribute::NONE;

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
	virtual void init(igl::opengl::glfw::Viewer *_viewer) override;
	void initialize();

	void launch();

	igl::opengl::ViewerData & mesh_by_id(int id);

	bool load();
	//bool load_param(std::string name);

	bool save();
	bool save_as(const std::string & save_as_dir);

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
	igl::ColorMapType cm = igl::ColorMapType::COLOR_MAP_TYPE_PARULA;
	double min_val = 0, max_val = 0;

	bool block_mouse_behavior(int button);
	void viewer_control();
	void viewer_control_2d();
	void viewer_control_3d();
	void draw_mesh();
	void fix_color(igl::opengl::ViewerData &data);
	void create_region_label();
	void build_region_edges(const Eigen::MatrixXd &pts, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2, Eigen::MatrixXd &C);

	///////////////
	// UI Panels //
	///////////////

	// Menu bar
	float draw_menu_bar(); // returns menu bar height

	// Left panel
	void draw_left_panel(float ypos, float width);
	void draw_file_menu();
	void draw_points_menu();
	void draw_mesh_menu();
	void draw_analysis_menu();

	// Right panel
	void draw_right_panel(float ypos, float width);
	void draw_histogram_menu();
	void draw_legend_menu();
	void draw_layer_menu();
	void draw_region_menu();

	// Toggle windows
	bool show_left_panel = true;
	bool show_file_menu = true;
	bool show_points_menu = true;
	bool show_mesh_menu = true;
	bool show_analysis_menu = true;

	bool show_right_panel = true;
	bool show_histogram_menu = true;
	bool show_legend_menu = true;
	bool show_layer_menu = true;
	bool show_region_menu = true;

public:
	// Menu stuff
	void draw_viewer_window() override;

	bool key_pressed(unsigned int unicode_key, int modifiers) override;
	bool key_up(int key, int modifiers) override;
};

// -----------------------------------------------------------------------------

} // namespace cellogram
