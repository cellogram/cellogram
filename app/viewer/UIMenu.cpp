////////////////////////////////////////////////////////////////////////////////
#include "FileDialog.h"
#include "UIState.h"
#include <cellogram/convex_hull.h>
#include <cellogram/laplace_energy.h>
#include <cellogram/mesh_solver.h>
#include <cellogram/region_grow.h>
#include <cellogram/remesh_adaptive.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex.h>
#include <cellogram/vertex_degree.h>
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

namespace {

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
		if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
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
		ImGui::InputFloat("hausd", &hausd);
		ShowTooltip("Maximal Hausdorff distance for the boundaries approximation.");
		ImGui::InputFloat("hgrad", &hgrad);
		ShowTooltip("Gradation value.");
		ImGui::InputFloat("hmax", &hmax);
		ShowTooltip("Maximal edge size.");
		ImGui::InputFloat("hmin", &hmin);
		ShowTooltip("Minimal edge size.");
		ImGui::InputFloat("hsiz", &hsiz);
		ShowTooltip("Build a constant size map of size x.");
		opt.hausd = hausd;
		opt.hgrad = hgrad;
		opt.hmax = hmax;
		opt.hmin = hmin;
		opt.hsiz = hsiz;
	}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

namespace {

	namespace AppLayout {
		constexpr float left_panel_width = 180;
		constexpr float right_panel_width = 180;
		constexpr float vertical_padding = 0;
		constexpr int height_colorbar = 20;
		constexpr int header_hue = 205;

		constexpr ImGuiWindowFlags window_flags =
			ImGuiWindowFlags_NoSavedSettings
			| ImGuiWindowFlags_AlwaysAutoResize;

		constexpr ImGuiTreeNodeFlags header_flags =
			ImGuiTreeNodeFlags_DefaultOpen
			| ImGuiTreeNodeFlags_OpenOnDoubleClick;
	};

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

void UIState::draw_viewer_window() {
	// Top menu bar
	float h = draw_menu_bar();

	// Menu on left
	if (show_left_panel) {
		draw_left_panel(h, menu_scaling() * AppLayout::left_panel_width);
	}

	// Menu on the right
	if (show_right_panel) {
		draw_right_panel(h, menu_scaling() * AppLayout::right_panel_width);
	}

	// Mess up with the mouse cursor
	if (delete_vertex || add_vertex) {
		// Cross hair
		ImGui::SetNextWindowPos(ImVec2(-100, -100), ImGuiSetCond_Always);
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
}

////////////////////////////////////////////////////////////////////////////////
// Main panels
////////////////////////////////////////////////////////////////////////////////

// Demonstrate creating a fullscreen menu bar and populating it.
float UIState::draw_menu_bar() {
	float h = 0;
	if (ImGui::BeginMainMenuBar()) {
		if (ImGui::BeginMenu("File")) {
			if (ImGui::MenuItem("Open", "Ctrl+O")) {
				std::string fname = FileDialog::openFileName(DATA_DIR, {"*.png", "*.tif", "*.tiff"});
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
			}

			if (ImGui::MenuItem("Quit", "Alt+F4")) {
				glfwSetWindowShouldClose(viewer.window, GLFW_TRUE);
			}
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("View")) {
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

			ImGui::EndMenu();
		}
		h = ImGui::GetWindowSize().y;
		ImGui::EndMainMenuBar();
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

	ImGui::SetNextWindowPos(ImVec2(0.0f, ypos), ImGuiSetCond_Always);
	ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_Always);
	ImGui::SetNextWindowSizeConstraints(ImVec2(width, 0.0f), ImVec2(width, height));

	ImGui::Begin("Pipeline", &show_left_panel, AppLayout::window_flags);

	if (SetHue(hue) && ImGui::CollapsingHeader("Input File", &show_file_menu, AppLayout::header_flags)) {
		draw_file_menu();
	}

	if (SetHue(hue) && ImGui::CollapsingHeader("Stage 1 - Detection", &show_points_menu, AppLayout::header_flags)) {
		draw_points_menu();
	}

	if (SetHue(hue) && ImGui::CollapsingHeader("Stage 2 - Matching", &show_mesh_menu, AppLayout::header_flags)) {
		draw_mesh_menu();
	}

	if (SetHue(hue) && ImGui::CollapsingHeader("Stage 3 - Analysis", &show_analysis_menu, AppLayout::header_flags)) {
		draw_analysis_menu();
	}

	ImGui::End();
}

// -----------------------------------------------------------------------------

void UIState::draw_right_panel(float ypos, float width) {
	const float hue = AppLayout::header_hue;

	auto canvas = ImGui::GetIO().DisplaySize;
	float xpos = canvas.x - width;
	float vpad = AppLayout::vertical_padding * menu_scaling();
	ypos += vpad;
	float height = canvas.y - ypos - vpad;

	ImGui::SetNextWindowPos(ImVec2(xpos, ypos), ImGuiSetCond_Always);
	ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiSetCond_FirstUseEver);
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

////////////////////////////////////////////////////////////////////////////////
// Left panel
////////////////////////////////////////////////////////////////////////////////

void UIState::draw_file_menu() {
	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;

	// Text filed showing loaded image
	std::vector<std::string> strings;
	std::istringstream f(save_dir);
	std::string s;
	while (std::getline(f, s, '\\')) {
		strings.push_back(s);
	}

	ImGui::Text("File: %s", strings.back().c_str());
	if (ImGui::IsItemHovered()) {
		ImGui::BeginTooltip();
		ImGui::PushTextWrapPos(w * 35.0f);
		ImGui::TextUnformatted(save_dir.c_str());
		ImGui::PopTextWrapPos();
		ImGui::EndTooltip();
	}

	if (ImGui::Button("Load Image", ImVec2(-1, 0))) {
		std::string fname = FileDialog::openFileName(DATA_DIR, {"*.png", "*.tif", "*.tiff"});
		if (!fname.empty()) {
			load_image(fname);

			show_image = true;

			// update UI
			viewer_control();
		}
	}

	if (state.mesh.points.size() == 0) {
		push_disabled();
	}
	if (ImGui::Button("Save##Points", ImVec2((w - p) / 2.f, 0))) {
		save();
	}
	if (state.mesh.points.size() == 0) {
		pop_disabled();
	}
    ImGui::SameLine(0, p);

	if (!data_available) {
		push_disabled();
	}
	if (ImGui::Button("Load##Points", ImVec2((w - p) / 2.f, 0))) {
		load();
	}
	if (!data_available) {
		pop_disabled();
	}
}

// -----------------------------------------------------------------------------

void UIState::draw_points_menu() {
	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
	ImGui::InputFloat("Sigma", &state.sigma, 0.1, 0, 2);
	ImGui::PopItemWidth();

	if (ImGui::Button("Detection", ImVec2((w - p), 0))) {
		detect_vertices();
	}

	bool was_delete = delete_vertex;
	if (was_delete) {
		push_selected();
	}
	if (ImGui::Button("Delete Vertex", ImVec2((w - p), 0))) {
		add_vertex = false;
		delete_vertex = !delete_vertex;
	}
	if (was_delete) {
		pop_selected();
	}

	bool was_add = add_vertex;
	if (was_add)
		push_selected();
	if (ImGui::Button("Add Vertex", ImVec2((w - p), 0))) {
		delete_vertex = false;
		add_vertex = !add_vertex;
	}
	if (was_add){
		pop_selected();
	}

	bool was_moved = move_vertex;
	if (was_moved)
		push_selected();
	if (ImGui::Button("Move vertex", ImVec2((w - p), 0))) {
		// drag a single vertex to a new starting position
		move_vertex = !move_vertex;
		viewer_control();
	}
	if (was_moved) {
		pop_selected();
	}
}

// -----------------------------------------------------------------------------

void UIState::draw_mesh_menu() {
	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);

	// disable if points are not detected
	if (state.mesh.points.size() == 0) {
		push_disabled();
	}
	ImGui::InputInt("Num Iter", &state.lloyd_iterations);

	if (ImGui::SliderFloat("t", &t, 0, 1)) {
		viewer_control();
	}

	ImGui::PopItemWidth();

	if (ImGui::Button("Untangle", ImVec2((w - p), 0))) {
		state.untangle();
		t = 1;
		mesh_color.resize(0, 0);
		viewer_control();
	}

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
	if (ImGui::Button("Lloyd", ImVec2((w - p), 0))) {
		state.relax_with_lloyd();
		if (!state.regions.empty()) {
			state.detect_bad_regions();
			state.check_regions();
			state.fix_regions();
		}
		t = 1;
		mesh_color.resize(0, 0);
		viewer_control();
	}
	// ImGui::PopItemWidth();

	ImGui::Checkbox("Fix", &state.fix_regular_regions);

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.50f);
	if (ImGui::SliderFloat("energy", &state.energy_variation_from_mean, 0, 5)) {
		selected_region = -1;
		state.detect_bad_regions();
		state.check_regions();

		create_region_label();

		viewer_control();
	}
	ImGui::PopItemWidth();

	if (ImGui::Button("build regions", ImVec2((w - p), 0))) {
		selected_region = -1;

		state.detect_bad_regions();
		state.check_regions();
		state.fix_regions();

		show_bad_regions = true;
		show_mesh_fill = true;

		create_region_label();

		viewer_control();
	}

	if (state.mesh.points.size() == 0) {
		pop_disabled();
	}

	// disable if regions are not availabe
	bool disable_region = state.regions.size() == 0;
	if (disable_region) {
		push_disabled();
	}

	if (ImGui::Button("solve regions", ImVec2((w - p), 0))) {
		state.resolve_regions();
		selected_region = -1;
		current_region_status = "";
		create_region_label();

		viewer_control();
	}
	if (disable_region) {
		pop_disabled();
	}

	if (ImGui::Button("Ultimate relaxation", ImVec2((w - p), 0))) {
		state.final_relax();
		t = 1;
		create_region_label();

		viewer_control();
	}
}

// -----------------------------------------------------------------------------

void UIState::draw_analysis_menu() {
	ImGui::Checkbox("Pillars", &state.image_from_pillars);

	auto reset_view_3d = [&]() {
		analysis_mode = true;
		show_mesh = false;
		show_image = false;
		show_mesh_fill = false;
		view_mode_3d = Mesh3DAttribute::NONE;
	};

	if (state.image_from_pillars) {
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.40f);
		ImGui::InputFloat("E [MPa]", &state.eps, 0.1, 0.01, 3);
		ShowTooltip("Young's modulus of pillars");
		ImGui::InputFloat("I [µm^4]", &state.I, 0.1, 0.01, 3);
		ShowTooltip("Area moment of inertia");
		ImGui::InputFloat("L [µm]", &state.L, 0.1, 0.01, 3);
		ShowTooltip("Length of pillars");
		ImGui::PopItemWidth(); //---> the resulting force is in micro-newtons (if displacements or in micrometers)

		if (ImGui::Button("Analyze 3D mesh", ImVec2(-1, 0))) {
			state.analyze_3d_mesh();
			reset_view_3d();
			view_mode_3d = Mesh3DAttribute::NORM_DISP;
			viewer_control();
		}
	} else {
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.40f);
		ImGui::InputFloat("Scaling [µm/px]", &state.scaling, 0.01, 0.001, 3);
		ImGui::InputFloat("Padding [µm]", &state.padding_size, 1, 0, 0);
		ImGui::InputFloat("Thickness [µm]", &state.thickness, 1, 0, 0);
		ImGui::InputFloat("Target volume (%)", &state.target_volume, 0, 0, 3);
		ShowTooltip("Target volume (for uniform 3D meshing only), in % of the bbox diagonal");
		ImGui::InputFloat("Power", &state.power, 0, 0, 3);
		ShowTooltip(
			"After the norm of the displacement field has been remapped to [0, 1]\n"
			"(0 being the largest displacement, 1 being no displacement), applies the\n"
			"power law x^p to the to the scalar field to produce the size map driving\nthe adaptive mesh.");
		ImGui::InputFloat2("Edge length (%)", state.target_mesh_size);
		ShowTooltip("Lower and upper bound of the size map driving the adaptive mesh");
		if (ImGui::TreeNode("Advanced meshing options")) {
			SetMmgOptions(state.mmg_options);
			ImGui::TreePop();
		}

		ImGui::InputFloat("E", &state.E, 0.1, 0.01, 3);
		ImGui::InputFloat("nu", &state.nu, 0.1, 0.01, 3);

		ImGui::PopItemWidth();

		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;

		if (ImGui::Button("Mesh 2D adaptive", ImVec2(2.5f*(w-p)/4.f, 0))) {
			state.mesh_2d_adaptive();
			reset_view_3d();
			viewer_control();
		}
		ImGui::SameLine(0, p);
		if (ImGui::Button("Extrude", ImVec2(1.5f*(w-p)/4.f, 0))) {
			state.extrude_2d_mesh();
			reset_view_3d();
			viewer_control();
		}
		ShowTooltip("Extrude the 2D mesh into a 3D volume mesh.");

		if (ImGui::Button("Mesh 3D uniform", ImVec2(2.5f*(w-p)/4.f, 0))) {
			state.mesh_3d_uniform();
			reset_view_3d();
			viewer_control();
		}
		ImGui::SameLine(0, p);
		if (ImGui::Button("Remesh", ImVec2(1.5f*(w-p)/4.f, 0))) {
			state.remesh_3d_adaptive();
			reset_view_3d();
			viewer_control();
		}
		ShowTooltip("Remesh the current 3D volume mesh adaptively using mmg3d.");

		ImGui::Spacing();
		ImGui::Separator();
		ImGui::Spacing();

		if (ImGui::Button("Analyze 3D mesh", ImVec2(-1, 0))) {
			state.analyze_3d_mesh();
			reset_view_3d();
			view_mode_3d = Mesh3DAttribute::NORM_DISP;
			viewer_control();
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Right panel
////////////////////////////////////////////////////////////////////////////////

void UIState::draw_histogram_menu() {
	if (hist.size() == 0) {
		compute_histogram();
	}

	const float hist_w = ImGui::GetWindowWidth() * 0.75f - 2;
	ImGui::PushItemWidth(hist_w + 2);

	static float min_img = 0;
	static float max_img = 1;

	auto pos = ImGui::GetWindowPos();
	int startX = pos.x + 10;
	int startY = pos.y + 52 * menu_scaling();

	ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
	ImGui::PlotHistogram("", hist.data(), hist.size(), 0, NULL, 0.0f, hist.maxCoeff(), ImVec2(0, 80));

	ImDrawList *draw_list = ImGui::GetWindowDrawList();
	draw_list->PushClipRectFullScreen();
	draw_list->AddLine(ImVec2(startX + hist_w * min_img, startY), ImVec2(startX + hist_w * min_img, startY + 75),
		IM_COL32(0, 255, 0, 255), 2.0f);
	draw_list->AddLine(ImVec2(startX + hist_w * max_img, startY), ImVec2(startX + hist_w * max_img, startY + 75),
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

		texture = (tmp.array() * 255).cast<unsigned char>();

		viewer_control();
	}
	if (ImGui::SliderFloat("max", &max_img, 0.f, 1.f)) {
		Eigen::MatrixXd tmp = (state.img.array() - min_img) / (max_img - min_img);
		tmp = tmp.unaryExpr(clamping);
		texture = (tmp.array() * 255).cast<unsigned char>();

		viewer_control();
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

	static GLuint color_bar_texture = -1;
	static const int width = ImGui::GetWindowWidth();
	if (color_bar_texture) {
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
		if (!analysis_mode) {
			viewer.core.trackball_angle = Eigen::Quaternionf::Identity();
		}
		viewer_control();
	}

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.30f);

	{
		auto labels = {"--", "X", "Y", "Z", "Norm"};
		auto tips = {"", "Displacement along X", "Displacement along Y", "Displacement along Z", "Norm of the displacement"};
		if (ComboWithTooltips("Show attribute", (int *)(&view_mode_3d), labels, tips)) {
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
