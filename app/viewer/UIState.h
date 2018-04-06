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
	int img_id;
	int bad_region_id;

	// UI options
	// double foo;

	bool image_loaded;


public:
	void initialize();

	void launch();

	igl::opengl::ViewerData & mesh_by_id(int id);

	virtual bool load(std::string name) override;

	virtual bool save(std::string name) override;

	void load_image(std::string name);
	void compute_hull();
	void compute_triangulation();

public:
	igl::opengl::ViewerData & points_data() { return mesh_by_id(points_id); }
	igl::opengl::ViewerData & hull_data() { return mesh_by_id(hull_id); }
	igl::opengl::ViewerData & img_data() { return mesh_by_id(img_id); }
	igl::opengl::ViewerData & bad_region_data() { return mesh_by_id(bad_region_id); }
private:
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
