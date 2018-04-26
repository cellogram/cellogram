#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "cellogram/State.h"
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

	int selected_region = -1;
	int selected_param = 0;
	// UI options
	// double foo;

	// Display flags
	float t;
	Eigen::RowVector4f vertex_color = Eigen::RowVector4f(1,0,0,0);
	Eigen::MatrixXd mesh_color;
	bool show_hull = false;
	//bool image_loaded = false;
	bool show_points = false;
	bool show_mesh_fill = true;
	bool show_image = false;
	bool show_matching = false;
	bool show_bad_regions = false;
	bool color_code = false;

	// Clicking flags
	bool select_region = false;
	bool add_vertex = false;
	bool delete_vertex = false;
	bool make_vertex_good = false;

	// Image
	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> img;

	std::string current_region_status;
public:
	void initialize();

	void launch();

	igl::opengl::ViewerData & mesh_by_id(int id);

	virtual bool load(std::string name) override;
	bool load_param(std::string name);

	virtual bool save(std::string name) override;

	virtual bool mouse_scroll(float delta_y) override;

	void load_image(std::string name);
	void detect_vertices();
	void display_image();
	void compute_hull();
	void clean_hull();
	void compute_triangulation();


	void reset_viewer();

public:
	igl::opengl::ViewerData & points_data() { return mesh_by_id(points_id); }
	igl::opengl::ViewerData & hull_data() { return mesh_by_id(hull_id); }
	igl::opengl::ViewerData & image_data() { return mesh_by_id(image_id); }
	igl::opengl::ViewerData & bad_region_data() { return mesh_by_id(bad_region_id); }
	igl::opengl::ViewerData & matching_data() { return mesh_by_id(matching_id); }
private:
	void viewer_control();
	void draw_mesh();
	void fix_color(igl::opengl::ViewerData &data);
	Eigen::VectorXd create_region_label();
	void build_region_edges(const Eigen::MatrixXd &pts, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2);

public:
	// Menu stuff
	void draw_viewer_menu() override;
	void draw_custom_window() override;
};

// -----------------------------------------------------------------------------

} // namespace cellogram
