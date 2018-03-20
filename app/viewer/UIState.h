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

	// UI options

public:
	void launch();

	virtual bool load(std::string name) override;

	virtual bool save(std::string name) override;

public:
	// Menu stuff
	void draw_viewer_menu() override;
	void draw_custom_window() override;
};

// -----------------------------------------------------------------------------

} // namespace cellogram
