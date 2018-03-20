#pragma once

// #include "State.hpp"

#include <igl/colormap.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

namespace cellogram {

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
	// State &state;

	// UI options

protected:
	void foo() {}
};

} // namespace cellogram

