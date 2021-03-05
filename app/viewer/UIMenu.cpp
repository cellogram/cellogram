////////////////////////////////////////////////////////////////////////////////
#include "FileDialog.h"
#include "UIState.h"

#include "IconsFontAwesome5.h"

#include <cellogram/convex_hull.h>
#include <cellogram/laplace_energy.h>
#include <cellogram/region_grow.h>
#include <cellogram/remesh_adaptive.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex.h>
#include <cellogram/vertex_degree.h>
#include <cellogram/image_reader.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Cylinder.h>

#include <igl/colormap.h>
#include <igl/jet.h>
#include <igl/opengl/gl.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/parula.h>
#include <igl/unproject_onto_mesh.h>

#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <GLFW/glfw3.h>
#include <algorithm>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

using zebrafish::logger;

	namespace {

		namespace AppLayout {
			constexpr float left_panel_width = 180;
			constexpr float right_panel_width = 180;
			constexpr float vertical_padding = 0;
			constexpr int height_colorbar = 20;
			constexpr int header_hue = 205;
			constexpr int button_height = 25;

			constexpr int arrow_button_size = 40;
			constexpr int icon_button_size = 28;

			static const float TIME_THRESHOLD = 0.5; //In secs

			constexpr ImGuiWindowFlags window_flags =
				ImGuiWindowFlags_NoSavedSettings
				| ImGuiWindowFlags_AlwaysAutoResize;

			constexpr ImGuiTreeNodeFlags header_flags =
				ImGuiTreeNodeFlags_DefaultOpen
				| ImGuiTreeNodeFlags_OpenOnDoubleClick;
		};


		void push_disabled() {
			ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
			ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
		}

		void pop_disabled() {
			ImGui::PopItemFlag();
			ImGui::PopStyleVar();
		}

		void push_selected() {
			ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyle().Colors[ImGuiCol_ButtonHovered]);
		}

		void pop_selected() {
			ImGui::PopStyleColor();
		}

		void push_color(float hue) {
			static float col_main_sat = 180.f/255.f;
			static float col_main_val = 161.f/255.f;
			static float col_area_sat = 124.f/255.f;
			static float col_area_val = 100.f/255.f;
			static float col_back_sat = 59.f/255.f;
			static float col_back_val = 40.f/255.f;
			ImVec4 col_text = ImColor::HSV(hue/255.f,  20.f/255.f, 235.f/255.f);
			ImVec4 col_main = ImColor::HSV(hue/255.f, col_main_sat, col_main_val);
			ImVec4 col_back = ImColor::HSV(hue/255.f, col_back_sat, col_back_val);
			ImVec4 col_area = ImColor::HSV(hue/255.f, col_area_sat, col_area_val);

			auto Vec4 = []( float r, float g, float b, float a ) {
				float h, s, v;
				ImGui::ColorConvertRGBtoHSV( r, g, b, h, s, v );
				ImGui::ColorConvertHSVtoRGB( h, s, v, r, g, b );
				return ImVec4(r,g,b,a);
			};

			ImGui::PushStyleColor(ImGuiCol_Header, Vec4(col_main.x, col_main.y, col_main.z, 0.56f));
			ImGui::PushStyleColor(ImGuiCol_HeaderHovered, Vec4(col_main.x, col_main.y, col_main.z, 0.86f));
			ImGui::PushStyleColor(ImGuiCol_HeaderActive, Vec4(col_main.x, col_main.y, col_main.z, 1.00f));
		};

		void pop_color() {
			ImGui::PopStyleColor();
			ImGui::PopStyleColor();
			ImGui::PopStyleColor();
		};

		struct SetHue {
			SetHue(float hue) { push_color(hue); }
			~SetHue() { pop_color(); }
			operator bool() const { return true; }
		};

		void draw_legend_item(float r, float g, float b, std::string label) {
			ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
			ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor(r / 255, g / 255, b / 255));
			ImGui::Button("    ");
			ImGui::SameLine();
			ImGui::Text("%s", label.c_str());
			ImGui::PopStyleColor(1);
			ImGui::PopItemFlag();
		}

		void ShowTooltip(const std::string &desc) {

			if ((ImGui::IsItemActive() || ImGui::IsItemHovered())  && GImGui->HoveredIdTimer > AppLayout::TIME_THRESHOLD) {
				ImGui::SetTooltip("%s", desc.c_str());
			}
		}

		bool ComboWithTooltips(const char *label, int *idx, std::vector<const char *> items,
			const std::vector<const char *> &tips) {
			bool value_changed = false;
			if (ImGui::BeginCombo(label, items[*idx])) {
				for (int i = 0; i < (int)items.size(); ++i) {
					bool is_selected = (*idx == i);
					if (ImGui::Selectable(items[i], is_selected)) {
						value_changed = true;
						*idx = i;
					}
					if (std::strlen(tips[i]) && (ImGui::IsItemActive() || ImGui::IsItemHovered())) {
						ImGui::SetTooltip("%s", tips[i]);
					}
					if (is_selected) {
						ImGui::SetItemDefaultFocus();
					}
				}
				ImGui::EndCombo();
			}
			return value_changed;
		}

		bool button_right(ImFont *icon_font)
		{
			ImGui::PushStyleVar(ImGuiStyleVar_ButtonTextAlign, ImVec2(0.5f, 0.6f));
			ImGui::PushFont(icon_font);
			bool ok = ImGui::Button(ICON_FA_CHEVRON_RIGHT, ImVec2(AppLayout::arrow_button_size, AppLayout::arrow_button_size));
			ImGui::PopFont();
			ImGui::PopStyleVar();
			ShowTooltip("Proceed to next step");
			return ok;
		}
		bool button_left(ImFont *icon_font)
		{
			ImGui::PushStyleVar(ImGuiStyleVar_ButtonTextAlign, ImVec2(0.5f, 0.6f));
			ImGui::PushFont(icon_font);
			bool ok = ImGui::Button(ICON_FA_CHEVRON_LEFT, ImVec2(AppLayout::arrow_button_size, AppLayout::arrow_button_size));
			ImGui::PopFont();
			ImGui::PopStyleVar();
			ShowTooltip("Go back to previous step");
			return ok;
		}
		bool icon_button(ImFont *icon_font, const char* button_name)
		{
			ImGui::PushStyleVar(ImGuiStyleVar_ButtonTextAlign, ImVec2(0.5f, 0.6f));
			ImGui::PushFont(icon_font);
			bool ok = ImGui::Button(button_name, ImVec2(AppLayout::icon_button_size, AppLayout::icon_button_size));
			ImGui::PopFont();
			ImGui::PopStyleVar();
			return ok;
		}

		void SetMmgOptions(MmgOptions &opt) {
		// -ar  x 	all codes 	Value for angle detection.
		// -hausd  x 	all codes 	Maximal Hausdorff distance for the boundaries approximation.
		// -hgrad  x 	all codes 	Gradation value.
		// -hmax  x 	all codes 	Maximal edge size.
		// -hmin  x 	all codes 	Minimal edge size.
		// -hsiz x 	all codes 	Build a constant size map of size x.
		// -noinsert 	all codes 	No point insertion/deletion.
		// -nomove 	all codes 	No point relocation.
		// -nosurf 	all codes 	No surface modifications.
		// -noswap 	all codes 	No edge flipping.
		// -nr 	all codes 	No angle detection.
		// -nsd  n 	mmg2d 	In mesh generation mode (no given triangle), save the subdomain of index n. Save all subdomains
		// if n=0 (default).

			float hausd = opt.hausd;
			float hgrad = opt.hgrad;
			float hmax = opt.hmax;
			float hmin = opt.hmin;
			float hsiz = opt.hsiz;
		// ImGui::InputFloat("hausd", &hausd);
		// ShowTooltip("Maximal Hausdorff distance for the boundaries approximation.");
			ImGui::InputFloat("hgrad", &hgrad);
			ShowTooltip("Gradation value.");
		// ImGui::InputFloat("hmax", &hmax);
		// ShowTooltip("Maximal edge size.");
		// ImGui::InputFloat("hmin", &hmin);
		// ShowTooltip("Minimal edge size.");
		// ImGui::InputFloat("hsiz", &hsiz);
		// ShowTooltip("Build a constant size map of size x.");
			opt.hausd = hausd;
			opt.hgrad = hgrad;
			opt.hmax = hmax;
			opt.hmin = hmin;
			opt.hsiz = hsiz;
		}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////


void UIState::draw_viewer_window() {

	// Top menu bar
	float h = draw_menu_bar();

	if(run_next == 5)
	{
		next_call_back();

		next_call_back = []() {};

		has_next_callback = false;
		run_next = 0;
		close_next = true;
	}

	if (has_next_callback)
	{
		run_next++;
	}


	// Menu on left
	if (show_left_panel) {
		draw_left_panel(h, menu_scaling() * AppLayout::left_panel_width);
	}

	// Menu on the right
	if (show_right_panel) {
		//draw_right_panel(h, menu_scaling() * AppLayout::right_panel_width);
	}

	// image viewer
	if (show_imageViewer_menu) {
		draw_image_viewer_menu();
	}

	// log window
	if (show_log_menu) {
		draw_log_window();
	}

    // property editor
    if (show_editor_menu) {
        draw_editor_window();
    }

	// Mess up with the mouse cursor
	if (delete_vertex || add_vertex) {
		// Cross hair
		ImGui::SetNextWindowPos(ImVec2(-100, -100), ImGuiCond_Always);
		ImGui::Begin("mouse_layer");
		ImVec2 p = ImGui::GetIO().MousePos;
		ImDrawList *draw_list = ImGui::GetWindowDrawList();
		draw_list->PushClipRectFullScreen();
		draw_list->AddLine(ImVec2(p.x - 50, p.y), ImVec2(p.x + 50, p.y),
			IM_COL32(delete_vertex ? 255 : 0, add_vertex ? 255 : 0, 0, 255), 2.0f);
		draw_list->AddLine(ImVec2(p.x, p.y - 50), ImVec2(p.x, p.y + 50),
			IM_COL32(delete_vertex ? 255 : 0, add_vertex ? 255 : 0, 0, 255), 2.0f);
		draw_list->PopClipRect();

		ImGui::GetIO().MouseDrawCursor = true;
		ImGui::SetMouseCursor(-1);
		ImGui::End();
	} else {
		ImGui::GetIO().MouseDrawCursor = false;
		ImGui::SetMouseCursor(0);
	}
	//ImGuiWindow* window = ImGui::GetCurrentWindow();
	//ImGuiContext& g = *GImGui;
	//
	//std::cout << "-----" << window->IDStack.size() << std::endl;
	//std::cout << window->DC.GroupStack.size() << std::endl;
	//std::cout << g.ColorModifiers.size() << std::endl;
	//std::cout << g.StyleModifiers.size() << std::endl;
	//std::cout << g.FontStack.size() << std::endl;

	UIsize.resize = false;  // after window resize, only redraw once
}

////////////////////////////////////////////////////////////////////////////////
// Main panels
////////////////////////////////////////////////////////////////////////////////

// Demonstrate creating a fullscreen menu bar and populating it.
float UIState::draw_menu_bar() {
	float h = 0;
	if (ImGui::BeginMainMenuBar()) {
		if (ImGui::BeginMenu("File")) {
			if (ImGui::MenuItem("Load image", "Ctrl+O")) {
				std::string fname = FileDialog::openFileName("./", {"*.tif", "*.tiff"});
				if (!fname.empty()) {
					load_image(fname);
					show_image = true;

					// update UI
					viewer_control();
				}
			}

			if (ImGui::MenuItem("Save", "Ctrl+S")) {
				save();
			}
			if (ImGui::MenuItem("Save As..")) {
				std::string fname = FileDialog::saveFileName("./", { "*.json" });
				if (!fname.empty()) {
					save_as(fname);
				}
			}

			if (ImGui::MenuItem("Quit", "Alt+F4")) {
				glfwSetWindowShouldClose(viewer.window, GLFW_TRUE);
			}
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("View")) {
			// Image
			ImGui::MenuItem("Image##Bar", nullptr, &show_image);

			ImGui::Separator();

			// Points
			ImGui::MenuItem("Points##Bar", nullptr, &show_points);
			ImGui::Indent();
			ImGui::ColorEdit3("Color", colorUI.mesh_vertex_color.data(),
				ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
			//ImGui::MenuItem("Coded##Bar", nullptr, &show_analysis_menu);

			ImGui::Unindent();
			ImGui::Separator();

			// Mesh
			ImGui::MenuItem("Mesh##Bar", nullptr, &show_mesh);
			ImGui::Indent();
			ImGui::ColorEdit3("Color", colorUI.mesh_line_color.data(),
				ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
			ImGui::MenuItem("Fill##Bar", nullptr, &show_mesh_fill);

			ImGui::Unindent();
			ImGui::Separator();

			// Analysis
			ImGui::MenuItem("Displacements##Bar", nullptr, false);
			ImGui::MenuItem("Tractions##Bar", nullptr, &show_traction_forces);

			ImGui::Separator();
			float current_pt_size = points_data().point_size;
			if(ImGui::SliderFloat("Point size", &current_pt_size, 0, 100))
			{
				points_data().point_size = current_pt_size;
				physical_data().point_size = current_pt_size;
			}

			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("Windows")) {
			/*
			// Left panel
			ImGui::MenuItem("Left Panel##Bar", nullptr, &show_left_panel);
			ImGui::Indent();
			ImGui::MenuItem("Input File##Bar", nullptr, &show_file_menu);
			ImGui::MenuItem("Detection##Bar", nullptr, &show_points_menu);
			ImGui::MenuItem("Matching##Bar", nullptr, &show_mesh_menu);
			ImGui::MenuItem("Analysis##Bar", nullptr, &show_analysis_menu);
			ImGui::Unindent();

			ImGui::Separator();

			// Right panel
			ImGui::MenuItem("Right Panel##Bar", nullptr, &show_right_panel);
			ImGui::Indent();
			ImGui::MenuItem("Histogram", nullptr, &show_histogram_menu);
			ImGui::MenuItem("Legend", nullptr, &show_legend_menu);
			ImGui::MenuItem("Layers", nullptr, &show_layer_menu);
			ImGui::MenuItem("Regions", nullptr, &show_region_menu);
			ImGui::Unindent();
			*/
			ImGui::MenuItem("Image Viewer", nullptr, &show_imageViewer_menu);
			ImGui::MenuItem("Log", nullptr, &show_log_menu);
            ImGui::MenuItem("Property Editor", nullptr, &show_editor_menu);

			ImGui::EndMenu();
		}
		// if (ImGui::BeginMenu("Regions")) {
		// 	// Left panel
		// 	ImGui::MenuItem("Select Region");
		// 	if (ImGui::MenuItem("Mark Good"))
		// 		std::cout << "bla bla" << std::endl;
		// 	ImGui::MenuItem("Mark bad");

		// 	ImGui::Separator();
		// 	ImGui::MenuItem("Split Region");
		// 	ImGui::MenuItem("Split Region");
		// 	ImGui::MenuItem("Split Region");

		// 	// Legend
		// 	ImGui::Separator();
		// 	ImGui::Text("Region Legend");
		// 	draw_legend_item(46, 204, 113, "Ok");
		// 	draw_legend_item(155, 89, 182, "Too Many Vertices");
		// 	draw_legend_item(241, 196, 15, "Too Few Vertices");
		// 	draw_legend_item(41, 128, 185, "Region Too Large");
		// 	draw_legend_item(192, 57, 43, "No Solution");
		// 	draw_legend_item(149, 165, 166, "Not Properly Closed");
		// 	draw_legend_item(52, 73, 94, "Other");

		// 	ImGui::EndMenu();
		// }

		if (false && ImGui::BeginMenu("Legend")) {
			draw_legend_item(46, 204, 113, "Ok");
			draw_legend_item(155, 89, 182, "Too Many Vertices");
			draw_legend_item(241, 196, 15, "Too Few Vertices");
			draw_legend_item(41, 128, 185, "Region Too Large");
			draw_legend_item(192, 57, 43, "No Solution");
			draw_legend_item(149, 165, 166, "Not Properly Closed");
			draw_legend_item(52, 73, 94, "Other");

			ImGui::EndMenu();
		}

		ImGui::SameLine(0, 100);
		ImGui::Text("%s", current_file_name.c_str());
		h = ImGui::GetWindowSize().y;
		UIsize.mainMenuHeight = ImGui::GetWindowHeight();
		ImGui::EndMainMenuBar();

		viewer_control();
	}
	return h;
}

// -----------------------------------------------------------------------------

void UIState::draw_left_panel(float ypos, float width) {
	const float hue = AppLayout::header_hue;

	auto canvas = ImGui::GetIO().DisplaySize;
	float vpad = AppLayout::vertical_padding * menu_scaling();
	ypos += vpad;
	float height = canvas.y - ypos - vpad;

	ImGui::SetNextWindowPos(ImVec2(0.0f, ypos), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiCond_Always);
	ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0.0f), ImVec2(width, height));

	if (true)
	{
		ImGui::Begin("Input File", NULL, AppLayout::window_flags);
		ypos += draw_file_menu();
		ImGui::End();
	}

	ypos += vpad;
	ImGui::SetNextWindowPos(ImVec2(0.0f, ypos), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiCond_Always);
	ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0.0f), ImVec2(width, height));

	if (state.phase_enumeration == 1)
	{
		ImGui::Begin("Stage 1 - Detection", &show_left_panel, AppLayout::window_flags);
		draw_points_menu();
		ImGui::End();
	}
	if (state.phase_enumeration == 2)
	{
		ImGui::Begin("Stage 2 - Meshing", &show_left_panel, AppLayout::window_flags);
		draw_mesh_menu();
		ImGui::End();
	}
	if (state.phase_enumeration == 3)
	{
		ImGui::Begin("Stage 3 - Depth", &show_left_panel, AppLayout::window_flags);
		draw_depth_menu();
		ImGui::End();
	}
	if (state.phase_enumeration == 4)
	{
		ImGui::Begin("Stage 4 - Analysis", &show_left_panel, AppLayout::window_flags);
		draw_analysis_menu();
		ImGui::End();
	}
	if (state.phase_enumeration == 5)
	{
		ImGui::Begin("Results", &show_left_panel, AppLayout::window_flags);
		draw_results_menu();
		ImGui::End();
	}
	//if (SetHue(hue) && state.phase_enumeration == 0 && ImGui::CollapsingHeader("Input File", &show_file_menu, AppLayout::header_flags)) {
	//	draw_file_menu();
	//}

	//if (SetHue(hue) && state.phase_enumeration == 1 && ImGui::CollapsingHeader("Detection", &show_points_menu, AppLayout::header_flags)) {
	//	draw_points_menu();
	//}

	//if (SetHue(hue) && state.phase_enumeration == 2 && ImGui::CollapsingHeader("Matching", &show_mesh_menu, AppLayout::header_flags)) {
	//	draw_mesh_menu();
	//}

	//if (SetHue(hue) && state.phase_enumeration == 3 && ImGui::CollapsingHeader("Analysis", &show_analysis_menu, AppLayout::header_flags)) {
	//	draw_analysis_menu();
	//}


}

// -----------------------------------------------------------------------------

void UIState::draw_right_panel(float ypos, float width) {
	const float hue = AppLayout::header_hue;

	auto canvas = ImGui::GetIO().DisplaySize;
	float xpos = canvas.x - width;
	float vpad = AppLayout::vertical_padding * menu_scaling();
	ypos += vpad;
	float height = canvas.y - ypos - vpad;

	ImGui::SetNextWindowPos(ImVec2(xpos, ypos), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiCond_FirstUseEver);
	ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0.0f), ImVec2(width, height));

	ImGui::Begin("Viewer", &show_right_panel, AppLayout::window_flags);

	if (SetHue(hue) && ImGui::CollapsingHeader("Histogram", &show_histogram_menu, AppLayout::header_flags)) {
		draw_histogram_menu();
	}

	if (SetHue(hue) && ImGui::CollapsingHeader("Legend", &show_legend_menu, AppLayout::header_flags)) {
		draw_legend_menu();
	}

	if (SetHue(hue) && ImGui::CollapsingHeader("Layers", &show_layer_menu, AppLayout::header_flags)) {
		draw_layer_menu();
	}

	if (SetHue(hue) && ImGui::CollapsingHeader("Region", &show_region_menu, AppLayout::header_flags)) {
		draw_region_menu();
	}

	ImGui::End();
}

// -----------------------------------------------------------------------------

void UIState::draw_image_viewer_menu() {

	if (UIsize.resize) {
		ImGui::SetNextWindowPos(ImVec2(UIsize.windowWidth - UIsize.rightWidth, UIsize.windowHeight - UIsize.imageViewerHeight));
		ImGui::SetNextWindowSize(ImVec2(UIsize.rightWidth, UIsize.imageViewerHeight));
	}

	if (!ImGui::Begin("Image Viewer", &show_imageViewer_menu)) {
		ImGui::End();
		return;
	}

	if (!current_file_name.empty() && !state.img3D.empty()) {
		// General Info
		ImGui::Text("Rows = %ld  Cols = %ld Layers = %lu", state.img3D[0].rows(), state.img3D[0].cols(), state.img3D.size());
		ImGui::PushItemWidth(UIsize.rightWidth / 2.0);
		// Select rotation type
		int rotation_type = static_cast<int>(viewer.core().rotation_type);
		static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
		static bool orthographic = true;
		if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0"))
		{
			using RT = igl::opengl::ViewerCore::RotationType;
			auto new_type = static_cast<RT>(rotation_type);
			if (new_type != viewer.core().rotation_type)
			{
				if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
				{
					trackball_angle = viewer.core().trackball_angle;
					orthographic = viewer.core().orthographic;
					viewer.core().trackball_angle = Eigen::Quaternionf::Identity();
					viewer.core().orthographic = true;
				}
				else if (viewer.core().rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
				{
					viewer.core().trackball_angle = trackball_angle;
					viewer.core().orthographic = orthographic;
				}
				viewer.core().set_rotation_type(new_type);
			}
		}

		if (ImGui::TreeNode("Advanced visualization")) {
			
			// Orthographic view
			ImGui::Checkbox("Orthographic projection", &(viewer.core().orthographic));

            // show marker radius
            if (ImGui::Checkbox("Show radius", &show_radiusPoints)) {
                if (show_radiusPoints) visual_data().point_size = 4;
                else visual_data().point_size = 6;  // shared by axis points
            }

			// z_multiplier
			ImGui::SliderFloat("Z-mult", &imgViewer.visual_z_mult, 1.0, 9.0, "%.1f");
			if (ImGui::IsItemHovered()) {
				ImGui::SetTooltip("Multiplier applied to z-axis for visualization in case the z-stack is too thin.");
			}

			ImGui::TreePop();
            ImGui::Separator();
        }

		ImGui::Separator(); ////////////////////////

		ImGui::PushItemWidth(UIsize.rightWidth / 2.0);
		std::vector<std::string> typeName{"Max Project", "Z-Slice"};
		ImGui::Combo("3D Image Viewer Type", &imgViewer.imageViewerType, typeName);
		ImGui::PopItemWidth();

		if (imgViewer.imageViewerType == 1) {
			ImGui::SliderInt("Slice", &imgViewer.sliceToShow, 0, state.img3D.size()-1);
		}

		ImGui::Spacing();

		if (ImGui::Checkbox("Show axis points", &show_axisPoints)) {
            visual_data().point_size = 6;
        }
		if (ImGui::IsItemHovered()) {
			ImGui::SetTooltip("Three array of points indicating the X, Y and Z axis");
		}
        ImGui::Checkbox("Show index", &show_allIndex);


		/*
		ImGui::Checkbox("Show axis points", &show_axisPoints);
		if (showTooltip && ImGui::IsItemHovered()) {
			ImGui::SetTooltip("Three array of points indicating the X, Y and Z axis");
		}
		*/
		ImGui::PopItemWidth();
	} else {
		ImGui::Text("No 3D image registered.");
	}

	ImGui::End();
}

// -----------------------------------------------------------------------------

void UIState::draw_log_window() {

	if (UIsize.resize) {
		ImGui::SetNextWindowPos(ImVec2(UIsize.leftWidth, UIsize.windowHeight - UIsize.bottomHeight));
		ImGui::SetNextWindowSize(ImVec2(UIsize.windowWidth - UIsize.leftWidth - UIsize.rightWidth, UIsize.bottomHeight));
	}

	if (!ImGui::Begin("Log", &show_log_menu))
	{
		ImGui::End();
		return;
	}

	ImGui::Separator();
	ImGui::BeginChild("scrolling", ImVec2(0, 0), false, ImGuiWindowFlags_HorizontalScrollbar);

	static ImGuiTextBuffer buf;

	// std::string log = oss.str();
	// ImGui::TextUnformatted(log.c_str());
	std::string log = oss.str();
	oss.str("");
	oss.clear();
	// AddLog(log.c_str(), buf);
	buf.appendf("%s", log.c_str());

	ImGui::TextUnformatted(buf.begin(), buf.end());

	if (ImGui::GetScrollY() >= ImGui::GetScrollMaxY())
		ImGui::SetScrollHereY(1.0f);

	ImGui::EndChild();
	ImGui::End();
}

// -----------------------------------------------------------------------------

// draw_editor_window() moved to UI_editor.cpp because it is too long

////////////////////////////////////////////////////////////////////////////////
// Left panel
////////////////////////////////////////////////////////////////////////////////

float UIState::draw_file_menu() {
	float h = 0;
	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;

	// Text filed showing loaded image
	// std::vector<std::string> strings;
	// std::istringstream f(save_dir);
	// std::string s;
	// while (std::getline(f, s, '\\')) {
	// 	strings.push_back(s);
	// }

	// ImGui::Text("File: %s", strings.back().c_str());
	if (ImGui::IsItemHovered()) {
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(w * 35.0f);
		ImGui::TextUnformatted(save_dir.c_str());
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.70f);
	if (ImGui::Button("Load Image", ImVec2(0, AppLayout::icon_button_size))) {
		std::string fname = FileDialog::openFileName("./", {"*.tif", "*.tiff" });
		if (!fname.empty()) {
			load_image(fname);
		}
	}

	ImGui::PopItemWidth();
	float spacing = ImGui::GetStyle().ItemInnerSpacing.x;
	ImGui::SameLine(0.0f, spacing);

	if (!data_available) {
		push_disabled();
	}
	if (icon_button(icon_font,ICON_FA_FOLDER_OPEN)) {
		load();
	}
	ShowTooltip("Load saved data for current image");
	if (!data_available) {
		pop_disabled();
	}

	ImGui::SameLine(0.0f, spacing);

	if (state.mesh.points.size() == 0) {
		push_disabled();
	}
	if (icon_button(icon_font, ICON_FA_SAVE)) {
		save();
	}
	ShowTooltip("Save data for current image");
	if (state.mesh.points.size() == 0) {
		pop_disabled();
	}
	h = ImGui::GetWindowSize().y;
	return h;
}

// -----------------------------------------------------------------------------
// Stage 1 menu

void UIState::draw_points_menu() {
	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
	ImGui::InputFloat("Sigma", &state.sigma, 0.1, 0, 2);
	ShowTooltip("Sigma for fitting Gaussians");
	ImGui::InputFloat("Otsu", &state.otsu_multiplier, 0.1, 0, 2);
	ShowTooltip("Multiplier to Otsu level for thresholding");
	ImGui::PopItemWidth();

	//if (ImGui::Button("Detection", ImVec2((w - p), 0))) {
	wait_window("wait_detect", "Detecting vertices", ICON_FA_COOKIE,
		[&]() {return ImGui::Button("Detect", ImVec2((w - p), 0));},
		[&]() {detect_vertices();});

		//viewer_control();
		//detect_vertices();
	//}
	// Histogram is only relevant for detection so should go here
	{
		ImGui::Separator();
		draw_histogram_menu();
	}


	if (ImGui::TreeNode("Advanced vertex options")) {
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.70f);
        ImGui::PushFont(icon_font);
		bool was_delete = delete_vertex;
		if (was_delete) {
			push_selected();
		}

		if (icon_button(icon_font, ICON_FA_TRASH_ALT)) {
			add_vertex = false;
			delete_vertex = !delete_vertex;
		}
		ShowTooltip("Delete points (d)");
		if (was_delete) {
			pop_selected();
		}

		bool was_add = add_vertex;
		if (was_add)
			push_selected();
		ImGui::SameLine();
		if (icon_button(icon_font, ICON_FA_PLUS)) {
			delete_vertex = false;
			add_vertex = !add_vertex;
		}
		ShowTooltip("Add missing points (a)");
		if (was_add) {
			pop_selected();
		}

		bool was_moved = move_vertex;
		if (was_moved)
			push_selected();
		ImGui::SameLine();
		if (icon_button(icon_font, ICON_FA_ARROWS_ALT)) {
			// drag a single vertex to a new starting position
			move_vertex = !move_vertex;
			viewer_control();
		}
		ShowTooltip("Move points to their correct position by clicking and dragging");
		if (was_moved) {
			pop_selected();
		}
        ImGui::PopFont();
		ImGui::PopItemWidth();
		ImGui::TreePop();
	}

	// Arrow buttons
	ImGui::Spacing();
	ImGui::Separator();
	ImGui::Spacing();
	{
		float spacing = ImGui::GetStyle().ItemInnerSpacing.x;
		if (button_left(icon_font))
		{
			add_vertex = false;
			move_vertex = false;
			delete_vertex = false;

			// add here also the clean up of this stage
			state.img.resize(0, 0);
			data_available = false;
			//phase_0();

			viewer_control();
		}
		ImGui::SameLine(0.0f, spacing);

		// disable if points are not detected
		const int nothing_detected = state.mesh.points.size();
		if (nothing_detected == 0) {
			push_disabled();
		}

		wait_window("wait_meshing", "Meshing", ICON_FA_PRAYING_HANDS,
			[&]() {return button_right(icon_font);},
			[&]()
		{
			add_vertex = false;
			move_vertex = false;
			delete_vertex = false;

			state.untangle();
			/*state.detect_bad_regions();
			state.check_regions();*/
			phase_2();
			viewer_control();
		});
		if (nothing_detected == 0) {
			pop_disabled();
		}
	}
}

// -----------------------------------------------------------------------------
// Stage 2 menu

void UIState::draw_mesh_menu() {
	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);

	// disable if points are not detected
	//if (state.mesh.points.size() == 0) {
	//	push_disabled();
	//}
	/*ImGui::InputInt("Num Iter", &state.lloyd_iterations);*/

	ImGui::PopItemWidth();

	//if (ImGui::Button("Untangle", ImVec2((w - p), 0))) {
	//	state.untangle();
	//	t = 1;
	//	mesh_color.resize(0, 0);
	//	viewer_control();
	//}

	//ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
	//if (ImGui::Button("Lloyd", ImVec2((w - p), 0))) {
	//	state.relax_with_lloyd();
	//	if (!state.regions.empty()) {
	//		state.detect_bad_regions();
	//		state.check_regions();
	//		// state.fix_regions();
	//	}
	//	t = 1;
	//	mesh_color.resize(0, 0);
	//	viewer_control();
	//}
	//ImGui::PopItemWidth();

	//ImGui::Checkbox("Fix", &state.fix_regular_regions);

	// ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
	// if (ImGui::SliderFloat("Energy", &state.energy_variation_from_mean, 0, 5)) {
	// 	selected_region = -1;
	// 	state.detect_bad_regions();
	// 	state.check_regions();

	// 	create_region_label();

	// 	viewer_control();
	// }
	// ImGui::PopItemWidth();
	// ShowTooltip("Set energy threshold for difficult to mesh regions");

	if(state.mesh.added_by_untangler.size() > 0) {
		if (ImGui::Button("Move vertex", ImVec2((w - p), 0))) {
		// drag a single vertex to a new starting position
			move_vertex = !move_vertex;
			viewer_control();
		}
	} else {
		ImGui::Text("Great. No vertex added.");
	}

	// if (ImGui::Button("Build regions", ImVec2((w - p), 0))) {
	// 	selected_region = -1;

	// 	state.detect_bad_regions();
	// 	state.check_regions();
	// 	// state.fix_regions();

	// 	show_bad_regions = true;
	// 	show_mesh_fill = true;

	// 	create_region_label();

	// 	viewer_control();
	// }
	// ShowTooltip("Detect difficult to mesh regions");

	//if (state.mesh.points.size() == 0) {
	//	pop_disabled();
	//}

	// disable if regions are not availabe
	// bool disable_region = state.regions.size() == 0;
	// if (disable_region) {
	// 	push_disabled();
	// }

	// //if (ImGui::Button("solve regions", ImVec2((w - p), 0))) {
	// //	state.resolve_regions();
	// //	selected_region = -1;
	// //	current_region_status = "";
	// //	create_region_label();

	// //	viewer_control();
	// //}

	// ShowTooltip("Solve difficult to mesh regions");
	// if (disable_region) {
	// 	pop_disabled();
	// }

	// if (ImGui::TreeNode("Advanced regions options"))
	// {
	// 	draw_region_menu();
	// 	ImGui::TreePop();
	// }

	// Arrow buttons
	ImGui::Spacing();
	ImGui::Separator();
	ImGui::Spacing();

	{
		float spacing = ImGui::GetStyle().ItemInnerSpacing.x;
		if (button_left(icon_font))
		{
			phase_1();
			// add here also the clean up of this stage
			state.mesh.triangles.resize(0, 0);
			state.mesh.added_by_untangler.resize(0);
			state.mesh.deleted_by_untangler.resize(0, 0);

			viewer_control();
		}
		ImGui::SameLine(0.0f, spacing);


		wait_window("wait_relax", "Relaxing mesh", ICON_FA_UMBRELLA_BEACH,
			[&]() {return button_right(icon_font);},
			[&]()
		{
			state.final_relax();
			phase_3();
			viewer_control();
		});
	}

}

// -----------------------------------------------------------------------------
// Stage 3 menu

void UIState::draw_depth_menu() {

const float width = ImGui::GetWindowWidth() * 0.5f;
ImGui::PushItemWidth(width);
    static double solverTol = 1e-8;
    if (ImGui::InputDouble("BSP Tol", &solverTol, 0.0, 0.0, "%.2e")) {
        state.bsp.Set_solverTol(solverTol);
    }

    if (ImGui::Button("Compute bspline")) {
        state.PrepareBsp();
    }
    if (ImGui::TreeNode("Advanced BSP")) {

        const int zN = state.img3D.size();
        const int xN = state.img3D[0].rows();
        const int yN = state.img3D[0].cols();
        static int xMin = 2;
        static int xMax = xN-3;
        static int xNum = xN-4;
        static int yMin = 2;
        static int yMax = yN-3;
        static int yNum = yN-4;
        static int zMin = 2;
        static int zMax = zN-3;
        static int zNum = zN-4;
        ImGui::SliderInt("xMin", &xMin, 2, xMax);
        ImGui::SliderInt("xMax", &xMax, xMin, xN-3);
        ImGui::InputInt("xNum", &xNum);
        ImGui::SliderInt("yMin", &yMin, 2, yMax);
        ImGui::SliderInt("yMax", &yMax, yMin, yN-3);
        ImGui::InputInt("yNum", &yNum);
        ImGui::SliderInt("zMin", &zMin, 2, zMax);
        ImGui::SliderInt("zMax", &zMax, zMin, zN-3);
        ImGui::InputInt("zNum", &zNum);
        if (ImGui::Button("Export interp result")) {  // DEBUG PURPOSE
            // new image
            std::vector<Eigen::MatrixXd> img;
            Eigen::VectorXd xArray = Eigen::VectorXd::LinSpaced(xNum, xMin, xMax);
            Eigen::VectorXd yArray = Eigen::VectorXd::LinSpaced(yNum, yMin, yMax);
            Eigen::VectorXd zArray = Eigen::VectorXd::LinSpaced(zNum, zMin, zMax);

            const auto InterpImage = [this, &xArray, &yArray, &zArray, &img]() {
                img.clear();
                for (int iz=0; iz<zArray.size(); iz++) {
                    Eigen::MatrixXd slice(xArray.size(), yArray.size());
                    for (int ix=0; ix<xArray.size(); ix++)
                        for (int iy=0; iy<yArray.size(); iy++) {
                            slice(ix, iy) = state.bsp.Interp3D(xArray(ix), yArray(iy), zArray(iz));
                        }
                    img.push_back(slice);
                }
            };
            InterpImage();
            bool succ = cellogram::WriteTif("/Users/ziyizhang/Projects/tmp/interp_cello.tif", img, 0, img.size()-1);
            logger().info("Saving result to interp_cello.tif. Status = {}", succ);
        }
        ShowTooltip("DEBUG PURPOSE. Will destroy underlying quadrature method.");
        ImGui::TreePop();
    }

    static float alpha = 0.5;
    if (ImGui::SliderFloat("Alpha", &alpha, 0.0, 1.0, "%.3f")) {
        zebrafish::cylinder::H = alpha;
    }
    static float defaultRadius = 4.0;
    if (ImGui::SliderFloat("Initial R", &defaultRadius, 2.0, 6.0)) {
        state.mesh.optimPara.defaultRadius = defaultRadius;
    }
    ImGui::Checkbox("Invert Color", &state.mesh.optimPara.invertColor);

    static int DSnum_round1 = state.img3D.size();
    ImGui::InputInt("DSnum round1", &DSnum_round1);
    static double DSgap_round1 = 0.5;
    ImGui::InputDouble("DSgap round1", &DSgap_round1);
    static double DSeps_round1 = 0.1;
    ImGui::InputDouble("DSeps round1", &DSeps_round1, 0.0, 0.0, "%.2e");
    static int DSnum_round2 = 10;
    ImGui::InputInt("DSnum round2", &DSnum_round2);
    static double DSgap_round2 = 0.1;
    ImGui::InputDouble("DSgap round2", &DSgap_round2);
    static double DSeps_round2 = 1e-4;
    ImGui::InputDouble("DSeps round2", &DSeps_round2, 0.0, 0.0, "%.2e");
    static bool updateInfo = false;
    ImGui::Checkbox("round 2 update info", &updateInfo);

ImGui::PopItemWidth();

    wait_window("wait_DC", "Depth Searching...", ICON_FA_UMBRELLA_BEACH,
        [&]() {return ImGui::Button("Depth Search", ImVec2(-1.0, 0));},
        [&]() {
        state.DepthSearch_FirstCall(DSnum_round1, DSgap_round1, DSeps_round1);
        state.DepthSearch_Refine(DSnum_round2, DSgap_round2, DSeps_round2, updateInfo);
    });

    if (ImGui::TreeNode("Advanced DC")) {
        if (ImGui::Button("Depth Search Refine")) {
            state.DepthSearch_Refine(DSnum_round2, DSgap_round2, DSeps_round2, updateInfo);
        }
        ImGui::TreePop();
    }

	if (ImGui::Button("To Next Stage"))
		phase_4();
}

// -----------------------------------------------------------------------------

void UIState::draw_analysis_menu() {

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.40f);
	if (ImGui::SliderFloat("Displacement", &imgViewer.deformScale, 0, 1)) {
		viewer_control();
	}
	ImGui::PopItemWidth();

	ImGui::Checkbox("Pillars", &state.image_from_pillars);

	ImGui::Spacing();
	ImGui::Separator();
	ImGui::Spacing();

	if (state.image_from_pillars) {
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
		ImGui::InputFloat("Scaling", &state.scaling, 0.01, 0.001, 3);
		ShowTooltip("Magnification factor of image [µm/px]");
		ImGui::InputFloat("E", &state.eps, 0.1, 0.01, 3);
		ShowTooltip("Young's modulus of pillars [MPa]");
		ImGui::InputFloat("I", &state.I, 0.1, 0.01, 3);
		ShowTooltip("Area moment of inertia [µm^4]");
		ImGui::InputFloat("L", &state.L, 0.1, 0.01, 3);
		ShowTooltip("Length of pillars [µm]");
		ImGui::PopItemWidth(); //---> the resulting force is in micro-newtons (if displacements or in micrometers)


		// Arrow buttons
		ImGui::Spacing();
		ImGui::Separator();
		ImGui::Spacing();

		{
			float spacing = ImGui::GetStyle().ItemInnerSpacing.x;
			if (button_left(icon_font))
			{
				phase_3();
				viewer_control();
			}

			ImGui::SameLine(0.0f, spacing);

			wait_window("wait_analysis_pillars", "Calculating bending forces", ICON_FA_COGS,
				[&]() {return button_right(icon_font);},
				[&]()
			{
				state.analyze_3d_mesh();
				analysis_mode = true;
				show_mesh = false;
				show_image = false;
				show_mesh_fill = false;

				phase_5();
				//reset_view_3d();
				view_mode_3d = Mesh3DAttribute::NORM_DISP;
				override_ranges = false;
				viewer_control();
			});
		}

	} else {

		// Material Model selection
		static int model_selection = 0;
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
		const char* material_models[] = { "Linear Elasticity", "Neo Hooke" };
		if (ImGui::Combo("Material Model", &model_selection, material_models, IM_ARRAYSIZE(material_models)))
		{
			switch (model_selection)
			{
			case 0: {
				state.formulation = "LinearElasticity";
				break;
			}
			case 1:
			{
				state.formulation = "NeoHookean";
				break;
			}
			default: std::cout << "invalid material model" << std::endl;
			}
		}

		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
		ImGui::InputFloat("Scaling", &state.scaling, 0.01, 0.001, 3);
		ShowTooltip("Magnification factor of image [µm/px]");
		ImGui::InputFloat("E", &state.E, 0.1, 0.01, 3);
		ShowTooltip("Young's modulus [kPa]");
		ImGui::InputFloat("nu", &state.nu, 0.1, 0.01, 3);
		ShowTooltip("Poisson's ratio");

		if (ImGui::TreeNode("Advanced mesh options"))
		{
			ImGui::InputFloat("Padding [µm]", &state.padding_size, 1, 0, 0);
			ImGui::InputFloat("Thickness [µm]", &state.thickness, 1, 0, 0);
			const char disp_abs[] = "Displ. Thres. [µm]";
			const char disp_rel[] = "Displ. Thres. [%]";
			ImGui::Checkbox("Relative Threshold", &state.relative_threshold);
			ShowTooltip("Set the threshold relative to the maximum displacement of the detected points.");
			ImGui::InputFloat(state.relative_threshold ? disp_rel : disp_abs, &state.displacement_threshold, 0, 0, 3);
			ShowTooltip("Threshold on the displacement to identify regions that will be meshed more finely.");
			ImGui::InputFloat("Edge Length", &state.uniform_mesh_size, 0, 0, 3);
			ShowTooltip(
				"Target edge length for the mesh generated for the physical simulation.\n"
				"If adaptive meshing is used, this value specifies the minimum edge length.\n"
				"The length is expressed in terms of the median edge-length of the\nreconstructed Delaunay mesh.");
			float hgrad = state.mmg_options.hgrad;
			ImGui::InputFloat("Gradation", &hgrad);
			ShowTooltip(
				"Use this parameter to control the ratio between the length of\n"
				"adjacent edges of the mesh generated for the physical simulation.");
			state.mmg_options.hgrad = std::max(hgrad, 1.f);
			ImGui::TreePop();
		}

		ImGui::PopItemWidth();

		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;

		//if (ImGui::Button("Mesh 2D adaptive", ImVec2(-1.0, 0))) {
		//	state.mesh_2d_adaptive();
		//	reset_view_3d();
		//	viewer_control();
		//	// Eigen::MatrixXd V;
		//	// Eigen::MatrixXi F;
		//	// Eigen::VectorXd S;
		//	// state.mesh.get_background_mesh(state.scaling, V, F, S, state.padding_size);
		//	// igl::colormap(cm, state.mesh.sizing, true, mesh_color);
		//	// V /= state.scaling;
		//	// points_data().clear();
		//	// points_data().show_faces = true;
		//	// points_data().set_mesh(V, F);
		//	// points_data().set_colors(mesh_color);
		//	// std::cout << state.mesh.sizing << std::endl;
		//	// physical_data().clear();
		//}

		//bool mesh_2d_empty = state.mesh3d.empty();
		//if (mesh_2d_empty) {
		//	push_disabled();
		//}
		//if (ImGui::Button("Mesh 3D volume", ImVec2(-1.0, 0))) {
		//	state.mesh_3d_volume();
		//	state.remesh_3d_adaptive();
		//	reset_view_3d();
		//	viewer_control();
		//}
		wait_window("wait_3d_mesh", "Generating 3D mesh for FE", ICON_FA_FIGHTER_JET,
			[&]() {return ImGui::Button("Build volumetric mesh", ImVec2(-1.0, 0));},
			[&]()
		{
			// [NOTE] zebrafish:
			// mesh_3d_volume():   simply call extrude_2d_mesh()
			// extrude_2d_mesh():  
			// 		- mesh_2d_adaptive():
			//				- get_background_mesh():     use "points" in Mesh to retrieve {mesh3d.V, mesh3d.F} and sizing scalar
			//				- remesh_adaptive_2d():		 use sizing scalar to re-mesh {mesh3d.V, mesh3d.F}
			//		- extrude_mesh():			 		 update the {mesh3d.V, mesh3d.F} from 3d surface to 3d tetrahedral mesh (inplace)
			// remesh_3d_adaptive():                     re-mesh {mesh3D.V, mesh3D.T} with sizing field

			state.mesh_3d_volume();
			state.remesh_3d_adaptive();

			analysis_mode = true;
			show_mesh = false;
			show_image = false;
			show_mesh_fill = false;
			view_mode_3d = Mesh3DAttribute::NONE;
			override_ranges = false;

			viewer_control();
		});

		// Arrow buttons
		ImGui::Spacing();
		ImGui::Separator();
		ImGui::Spacing();

		{
			float spacing = ImGui::GetStyle().ItemInnerSpacing.x;
			if (button_left(icon_font))
			{
				phase_3();
				viewer_control();
			}

			bool mesh_3d_empty = state.mesh3d.empty();
			if (mesh_3d_empty) {
				push_disabled();
			}

			ImGui::SameLine(0.0f, spacing);

			//if (button_right(icon_font))
			//{
			//	state.analyze_3d_mesh();
			//	//state.phase_enumeration = 4;
			//	phase_4();
			//	//reset_view_3d();
			//	view_mode_3d = Mesh3DAttribute::NORM_DISP;
			//	viewer_control();
			//}
			wait_window("wait_analysis", "Running FEA", ICON_FA_COGS,
				[&]() {return button_right(icon_font);},
				[&]()
			{
				state.analyze_3d_mesh();
				phase_5();
				//reset_view_3d();
				view_mode_3d = Mesh3DAttribute::NORM_DISP;
				override_ranges = false;
				viewer_control();

				// TODO: add also the pillar calculation
			});
			if (mesh_3d_empty) {
				pop_disabled();
			}
		}

	}
}

void UIState::draw_results_menu()
{
	// if(!state.image_from_pillars)
	{
		static int view_current = 0;
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.4f);
		// const char* views[] = { "U", "S", "T"};
		const char* views[] = { "U", "T"};
		ImGui::Combo("View", &view_current, views, IM_ARRAYSIZE(views));

		static int sub_view_current = 0;
		if (view_current == 0){
			if(show_traction_forces)
				override_ranges = false;
			show_traction_forces = false;
			const char* sub_views[] = { "Mag", "Ux", "Uy", "Uz" };
			if (ImGui::Combo("Uview##View", &sub_view_current, sub_views, IM_ARRAYSIZE(sub_views)))
			{
				Mesh3DAttribute new_view = Mesh3DAttribute::NONE;

				switch (sub_view_current)
				{
				case 0: new_view = Mesh3DAttribute::NORM_DISP; break;
				case 1: new_view = Mesh3DAttribute::X_DISP; break;
				case 2: new_view = Mesh3DAttribute::Y_DISP; break;
				case 3: new_view = Mesh3DAttribute::Z_DISP; break;
				default:
					assert(false);
				}
				if(view_mode_3d != new_view)
					override_ranges = false;
				view_mode_3d = new_view;
			};
		}
		// else if (view_current == 1){
		// 	const char* sub_views[] = { "Mises", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz" };
		// 	ImGui::Combo("Sview##View", &sub_view_current, sub_views, IM_ARRAYSIZE(sub_views));
		// }
		else if (view_current == 1){
			if (!show_traction_forces)
				override_ranges = false;
			show_traction_forces = true;
			const char* sub_views[] = { "Mag", "Tx", "Ty", "Tz" };
			if(ImGui::Combo("Tview##View", &sub_view_current, sub_views, IM_ARRAYSIZE(sub_views)));
			{
				Mesh3DAttribute new_view = Mesh3DAttribute::NONE;

				switch (sub_view_current)
				{
				case 0: new_view = Mesh3DAttribute::NORM_DISP; break;
				case 1: new_view = Mesh3DAttribute::X_DISP; break;
				case 2: new_view = Mesh3DAttribute::Y_DISP; break;
				case 3: new_view = Mesh3DAttribute::Z_DISP; break;
				default:
					assert(false);
				}

				if (view_mode_3d != new_view)
					override_ranges = false;
				view_mode_3d = new_view;
			}
		}
		ImGui::PopItemWidth();
		if(!state.image_from_pillars && show_traction_forces)
			ImGui::Checkbox("Smooth results", &show_smoothed_results);
	}

	// Colorbar
	static GLuint color_bar_texture = 0;
	const int width = ImGui::GetWindowWidth();
	if (color_bar_texture == 0) {
		Eigen::Matrix<unsigned char, Eigen::Dynamic, 4, Eigen::RowMajor> cmap(width * AppLayout::height_colorbar, 4);

		Eigen::MatrixXd t = Eigen::VectorXd::LinSpaced(width, 0, width);
		Eigen::MatrixXd col;
		igl::colormap(cm, t, true, col);
		assert(col.rows() == width);
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < AppLayout::height_colorbar; ++j) {
				for (int c = 0; c < 3; ++c)
					cmap(j * width + i, c) = col(i, c) * 255;
			}
		}
		cmap.col(3).setConstant(255);

		glGenTextures(1, &color_bar_texture);
		glBindTexture(GL_TEXTURE_2D, color_bar_texture);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, AppLayout::height_colorbar, 0, GL_RGBA, GL_UNSIGNED_BYTE, cmap.data());
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	ImGui::Image(reinterpret_cast<ImTextureID>(color_bar_texture), ImVec2(width, AppLayout::height_colorbar));

	if (std::abs(min_val) <= 1e-20)
		ImGui::Text("0");
	else {
		const int min_power = floor(log10(std::abs(min_val)));
		ImGui::Text("%ge%d", round(min_val * pow(10, -min_power) * 100) / 100., min_power);
	}

	if (std::abs(max_val) <= 1e-20) {
		ImGui::SameLine(width - 10);
		ImGui::Text("0");
	}
	else {
		ImGui::SameLine(width - 40);
		const int max_power = floor(log10(std::abs(max_val)));
		ImGui::Text("%ge%d", round(max_val * pow(10, -max_power) * 100) / 100., max_power);
	}

	ImGui::PushItemWidth(ImGui::GetWindowWidth()*0.75);

	float tmpmin = min_val;
	if (ImGui::SliderFloat("min##Color", &tmpmin, real_min_val, real_max_val))
	{
		override_ranges=true;
		min_val = tmpmin;
		viewer_control();
	}
	float tmpmax = max_val;
	if (ImGui::SliderFloat("max##Color", &tmpmax, real_min_val, real_max_val))
	{
		override_ranges = true;
		max_val = tmpmax;
		viewer_control();
	}
	ImGui::PopItemWidth();

	// Arrow buttons
	ImGui::Spacing();
	ImGui::Separator();
	ImGui::Spacing();

	float spacing = ImGui::GetStyle().ItemInnerSpacing.x;
	ImGui::PushFont(icon_font);
	if (button_left(icon_font))
	{
		state.mesh3d.clear();
		// go back to analysis stage and possibly remove current solution
		phase_4();
		viewer_control();
	}
	ImGui::SameLine(0.0f, spacing);
	if (ImGui::Button(ICON_FA_SAVE, ImVec2(AppLayout::arrow_button_size, AppLayout::arrow_button_size)))
	{
		// save solution
		save();
	}
	ImGui::PopFont();
}


// -----------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////
// Right panel
////////////////////////////////////////////////////////////////////////////////

void UIState::draw_histogram_menu() {
	if (hist.size() == 0) {
		compute_histogram();
	}

	const float hist_w = ImGui::GetWindowWidth() * 0.75f - 2;
	const float hist_h = 80 * menu_scaling();
	ImGui::PushItemWidth(hist_w + 2);

	static float min_img = 0;
	static float max_img = 1;

	auto before = ImGui::GetCursorScreenPos();
	ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
	ImGui::PlotHistogram("", hist.data(), hist.size(), 0, NULL, 0.0f, hist.maxCoeff(), ImVec2(0, hist_h));
	auto after = ImGui::GetCursorScreenPos();
	after.y -= ImGui::GetStyle().ItemSpacing.y;

	ImDrawList *draw_list = ImGui::GetWindowDrawList();
	draw_list->PushClipRectFullScreen();
	draw_list->AddLine(ImVec2(before.x + hist_w * min_img, before.y),
		ImVec2(before.x + hist_w * min_img, after.y),
		IM_COL32(0, 255, 0, 255), 2.0f);
	draw_list->AddLine(ImVec2(before.x + hist_w * max_img, before.y),
		ImVec2(before.x + hist_w * max_img, after.y),
		IM_COL32(0, 255, 0, 255), 2.0f);
	draw_list->PopClipRect();

	ImGui::PopItemFlag();

	static const auto clamping = [](const double x) {
		if (x > 0) {
			if (x < 1)
				return x;
			else
				return 1.;
		} else {
			return 0.0;
		}
	};

	if (ImGui::SliderFloat("min", &min_img, 0.f, 1.f)) {
		Eigen::MatrixXd tmp = (state.img.array() - min_img) / (max_img - min_img);
		tmp = tmp.unaryExpr(clamping);
		texture = (tmp.array() * 255 * imgViewer.darkenScale).cast<unsigned char>();

		viewer_control();
	}
	if (ImGui::SliderFloat("max", &max_img, 0.f, 1.f)) {
		Eigen::MatrixXd tmp = (state.img.array() - min_img) / (max_img - min_img);
		tmp = tmp.unaryExpr(clamping);
		texture = (tmp.array() * 255 * imgViewer.darkenScale).cast<unsigned char>();

		viewer_control();
	}
	if (ImGui::SliderFloat("dark", &imgViewer.darkenScale, 0.01, 1.0)) {
		Eigen::MatrixXd tmp = (state.img.array() - min_img) / (max_img - min_img);
		tmp = tmp.unaryExpr(clamping);
		texture = (tmp.array() * 255 * imgViewer.darkenScale).cast<unsigned char>();

		viewer_control();
	}
	if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Reduce the intensity of bright colors.\nNo effect on results, only for visualization purpose.");
    }

	ImGui::PopItemWidth();
}

// -----------------------------------------------------------------------------

void UIState::draw_legend_menu() {
	draw_legend_item(46, 204, 113, "Ok");
	draw_legend_item(155, 89, 182, "Too Many Vertices");
	draw_legend_item(241, 196, 15, "Too Few Vertices");
	draw_legend_item(41, 128, 185, "Region Too Large");
	draw_legend_item(192, 57, 43, "No Solution");
	draw_legend_item(149, 165, 166, "Not Properly Closed");
	draw_legend_item(52, 73, 94, "Other");

	ImGui::Separator();

	static GLuint color_bar_texture = 0;
	static const int width = ImGui::GetWindowWidth();
	if (color_bar_texture == 0) {
		Eigen::Matrix<unsigned char, Eigen::Dynamic, 4, Eigen::RowMajor> cmap(width * AppLayout::height_colorbar, 4);

		Eigen::MatrixXd t = Eigen::VectorXd::LinSpaced(width, 0, width);
		Eigen::MatrixXd col;
		igl::colormap(cm, t, true, col);
		assert(col.rows() == width);
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < AppLayout::height_colorbar; ++j) {
				for (int c = 0; c < 3; ++c)
					cmap(j * width + i, c) = col(i, c) * 255;
			}
		}
		cmap.col(3).setConstant(255);

		glGenTextures(1, &color_bar_texture);
		glBindTexture(GL_TEXTURE_2D, color_bar_texture);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, AppLayout::height_colorbar, 0, GL_RGBA, GL_UNSIGNED_BYTE, cmap.data());
		glGenerateMipmap(GL_TEXTURE_2D);
		glBindTexture(GL_TEXTURE_2D, 0);
	}
	ImGui::Image(reinterpret_cast<ImTextureID>(color_bar_texture), ImVec2(width, AppLayout::height_colorbar));

	if (std::abs(min_val) <= 1e-20)
		ImGui::Text("0");
	else {
		const int min_power = floor(log10(std::abs(min_val)));
		ImGui::Text("%ge%d", round(min_val * pow(10, -min_power) * 100) / 100., min_power);
	}

	if (std::abs(max_val) <= 1e-20) {
		ImGui::SameLine(width - 10);
		ImGui::Text("0");
	} else {
		ImGui::SameLine(width - 40);
		const int max_power = floor(log10(std::abs(max_val)));
		ImGui::Text("%ge%d", round(max_val * pow(10, -max_power) * 100) / 100., max_power);
	}
}

// -----------------------------------------------------------------------------

void UIState::draw_layer_menu() {
	if (ImGui::Checkbox("", &show_mesh)) {
		viewer_control();
	}
	ImGui::SameLine();
	if (ImGui::ColorEdit4("Mesh", points_data().line_color.data(),
		ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel)) {
		viewer_control();
}

if (ImGui::Checkbox("Show hull", &show_hull)) {
	viewer_control();
}
if (ImGui::Checkbox("Points", &show_points)) {
	viewer_control();
}
ImGui::SameLine();
if (ImGui::Checkbox("Coded", &color_code)) {
	viewer_control();
}

if (ImGui::Checkbox("Mesh Fill", &show_mesh_fill)) {
	viewer_control();
}
if (ImGui::Checkbox("Show image", &show_image)) {
	viewer_control();
}
if (ImGui::Checkbox("Show matching", &show_matching)) {
	viewer_control();
}
if (ImGui::Checkbox("Bad regions", &show_bad_regions)) {
	viewer_control();
}
if (ImGui::Checkbox("Selected region", &show_selected_region)) {
	viewer_control();
		//debugging region, can be deleted after
	state.check_region(selected_region);
}

if (ImGui::Checkbox("Enable 3D view", &analysis_mode)) {
	/*
	if (!analysis_mode) {
		viewer.core().trackball_angle = Eigen::Quaternionf::Identity();
	}
	*/
	viewer_control();
}

ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.30f);

{
	auto labels = {"--", "X", "Y", "Z", "Norm"};
	auto tips = {"", "Displacement along X", "Displacement along Y", "Displacement along Z", "Norm of the displacement"};
	if (ComboWithTooltips("Show attribute", (int *)(&view_mode_3d), labels, tips)) {
		override_ranges = false;
		viewer_control();
	}
}

if (ImGui::Checkbox("Traction forces", &show_traction_forces)) {
	viewer_control();
}

ImGui::PopItemWidth();
}

////////////////////////////////////////////////////////////////////////////////
// Toolbox
////////////////////////////////////////////////////////////////////////////////

void UIState::draw_region_menu() {
	if (state.regions.size() == 0) {
		push_disabled();
	}
	bool region_was_selected = select_region;
	if (region_was_selected) {
		push_selected();
	}
	if (ImGui::Button("Select Region", ImVec2(-1, 0))) {
		deselect_all_buttons();
		create_region_label();
		select_region = true;
		show_selected_region = true;
	}
	if (region_was_selected) {
		pop_selected();
	}
	if (state.regions.size() == 0) {
		pop_disabled();
	}

	if (state.regions.size() == 0) {
		push_disabled();
	}
	if (ImGui::Button("Export Region", ImVec2(-1, 0))) {
		deselect_all_buttons();
		export_region();
	}
	if (state.regions.size() == 0) {
		pop_disabled();
	}

	int n_was_selected = selected_region;
	if (n_was_selected < 0) {
		push_disabled();
	}
	if (ImGui::Button("Grow Selected", ImVec2(-1, 0))) {
		state.grow_region(selected_region);
		viewer_control();
	}

	if (ImGui::Button("Solve Selected", ImVec2(-1, 0))) {
		state.resolve_region(selected_region);

		selected_region = -1;
		current_region_status = "";
		viewer_control();
	}
	if (n_was_selected < 0) {
		pop_disabled();
	}

	ImGui::Separator();

	bool was_marked_good = make_vertex_good;
	if (was_marked_good) {
		push_selected();
	}
	if (ImGui::Button("Mark good", ImVec2(-1, 0))) {
		deselect_all_buttons();
		// select vertices and mark them as good permanently
		make_vertex_good = true;
		viewer_control();
	}
	if (was_marked_good) {
		pop_selected();
	}

	bool was_marked_bad = make_vertex_bad;
	if (was_marked_bad) {
		push_selected();
	}
	if (ImGui::Button("Mark bad", ImVec2(-1, 0))) {
		deselect_all_buttons();
		// select vertices and mark them as good permanently
		make_vertex_bad = true;
		viewer_control();
	}
	if (was_marked_bad) {
		pop_selected();
	}

	ImGui::Separator();

	if (state.regions.size() == 0) {
		push_disabled();
	}
	bool was_split = split_region > -1;
	if (was_split) {
		push_selected();
	}
	if (ImGui::Button("Split Region", ImVec2(-1, 0))) {
		split_region = 0;
		deselect_all_buttons();
		// split_region = true;
		viewer_control();
	}
	if (was_split) {
		pop_selected();
	}
	if (state.regions.size() == 0) {
		pop_disabled();
	}

	if (selected_region >= 0) {
		ImGui::Separator();
		ImGui::PushItemWidth(ImGui::GetWindowWidth());
		int nVi = state.regions[selected_region].region_interior.size();
		int nVtotal = state.regions[selected_region].region_boundary.size() + nVi;
		int nTri = state.regions[selected_region].region_triangles.size();
		ImGui::LabelText("", "Region %d ", selected_region);
		if (state.regions[selected_region].points_delta != 0)
			ImGui::LabelText("", "%s (%i) ", current_region_status.c_str(),
				state.regions[selected_region].points_delta);
		else
			ImGui::LabelText("", "%s ", current_region_status.c_str());
		ImGui::LabelText("", "#V: %i (%i) ", nVtotal, nVi);
		ImGui::LabelText("", "#F: %i ", nTri);
		ImGui::PopItemWidth();
	}
}

// -----------------------------------------------------------------------------

} // namespace cellogram
