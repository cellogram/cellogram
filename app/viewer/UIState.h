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

	// UI options
	// double foo;

public:
	void initialize();

	void launch();

	igl::opengl::ViewerData & mesh_by_id(int id);

	virtual bool load(std::string name) override;

	virtual bool save(std::string name) override;

	void compute_hull();

public:
	igl::opengl::ViewerData & points_data() { return mesh_by_id(points_id); }
	igl::opengl::ViewerData & hull_data() { return mesh_by_id(hull_id); }

public:
	// Menu stuff
	void draw_viewer_menu() override;
	void draw_custom_window() override;
};

// -----------------------------------------------------------------------------

} // namespace cellogram
