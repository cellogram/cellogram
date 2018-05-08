////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include "FileDialog.h"
#include <cellogram/laplace_energy.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex.h>
#include <cellogram/region_grow.h>
#include <cellogram/vertex_degree.h>
#include <cellogram/mesh_solver.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <imgui/imgui_internal.h>
#include <imgui/imgui.h>
#include <algorithm>
#include <vector>
#include <igl/parula.h>
#include <igl/jet.h>
#include <cellogram/convex_hull.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

	namespace
	{
		void push_disabled()
		{
			ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
			ImGui::PushStyleVar(ImGuiStyleVar_Alpha, ImGui::GetStyle().Alpha * 0.5f);
		}

		void pop_disabled()
		{
			ImGui::PopItemFlag();
			ImGui::PopStyleVar();
		}

		void push_selected()
		{
			ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyle().Colors[ImGuiCol_ButtonHovered]);
		}

		void pop_selected()
		{
			ImGui::PopStyleColor();
		}

		void draw_legend_item(float r, float g, float b, std::string label)
		{
			ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
			ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor(r / 255, g / 255, b / 255));
			ImGui::Button("    ");
			ImGui::SameLine();
			ImGui::Text(label.c_str());
			ImGui::PopStyleColor(1);
			ImGui::PopItemFlag();
		}
	}

// -----------------------------------------------------------------------------

void UIState::draw_viewer_menu() {
	//// Viewing options
	//if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen)) {
	//	if (ImGui::Button("Center object", ImVec2(-1, 0))) {
	//		viewer.core.align_camera_center(viewer.data().V, viewer.data().F);
	//	}
	//	if (ImGui::Button("Snap canonical view", ImVec2(-1, 0))) {
	//		viewer.snap_to_canonical_quaternion();
	//	}

	//	// Zoom
	//	ImGui::PushItemWidth(80 * menu_scaling());
	//	ImGui::DragFloat("Zoom", &(viewer.core.camera_zoom), 0.05f, 0.1f, 20.0f);

	//	// Select rotation type
	//	static int rotation_type = static_cast<int>(viewer.core.rotation_type);
	//	static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
	//	static bool orthographic = true;
	//	if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0")) {
	//		using RT = igl::opengl::ViewerCore::RotationType;
	//		auto new_type = static_cast<RT>(rotation_type);
	//		if (new_type != viewer.core.rotation_type) {
	//			if (new_type == RT::ROTATION_TYPE_NO_ROTATION) {
	//				trackball_angle = viewer.core.trackball_angle;
	//				orthographic = viewer.core.orthographic;
	//				viewer.core.trackball_angle = Eigen::Quaternionf::Identity();
	//				viewer.core.orthographic = true;
	//			} else if (viewer.core.rotation_type == RT::ROTATION_TYPE_NO_ROTATION) {
	//				viewer.core.trackball_angle = trackball_angle;
	//				viewer.core.orthographic = orthographic;
	//			}
	//			viewer.core.set_rotation_type(new_type);
	//		}
	//	}

	//	// Orthographic view
	//	ImGui::Checkbox("Orthographic view", &(viewer.core.orthographic));
	//	ImGui::PopItemWidth();
	//}

	// Draw options
	//if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen)) {
	//	if (ImGui::Checkbox("Face-based", &(viewer.data().face_based))) {
	//		viewer.data().set_face_based(viewer.data().face_based);
	//	}
	//	ImGui::Checkbox("Show texture", &(viewer.data().show_texture));
	//	if (ImGui::Checkbox("Invert normals", &(viewer.data().invert_normals))) {
	//		viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
	//	}
	//	ImGui::Checkbox("Show overlay", &(viewer.data().show_overlay));
	//	ImGui::Checkbox("Show overlay depth", &(viewer.data().show_overlay_depth));
	//	ImGui::ColorEdit4("Background", viewer.core.background_color.data(),
	//			ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
	//	ImGui::ColorEdit4("Line color", viewer.data().line_color.data(),
	//			ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
	//	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
	//	ImGui::DragFloat("Shininess", &(viewer.data().shininess), 0.05f, 0.0f, 100.0f);
	//	ImGui::PopItemWidth();
	//}

	//// Overlays
	//if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen)) {
	//	ImGui::Checkbox("Wireframe", &(viewer.data().show_lines));
	//	ImGui::Checkbox("Fill", &(viewer.data().show_faces));
	//	ImGui::Checkbox("Show vertex labels", &(viewer.data().show_vertid));
	//	ImGui::Checkbox("Show faces labels", &(viewer.data().show_faceid));
	//}
}

// Demonstrate creating a fullscreen menu bar and populating it.
void UIState::build_menu_bar()
{
	if (ImGui::BeginMainMenuBar())
	{
		if (ImGui::BeginMenu("File"))
		{
			if (ImGui::MenuItem("Open", "Alt+O")) {
				std::string fname = FileDialog::openFileName(DATA_DIR, { "*.png", "*.tif", "*.tiff" });
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
			if (ImGui::MenuItem("!!Save As..")) {}
			ImGui::Separator();
			if (ImGui::BeginMenu("Options"))
			{
				static bool enabled = true;
				ImGui::MenuItem("Enabled", "", &enabled);
				ImGui::BeginChild("child", ImVec2(0, 60), true);
				for (int i = 0; i < 10; i++)
					ImGui::Text("Scrolling Text %d", i);
				ImGui::EndChild();
				static float f = 0.5f;
				static int n = 0;
				static bool b = true;
				ImGui::SliderFloat("Value", &f, 0.0f, 1.0f);
				ImGui::InputFloat("Input", &f, 0.1f);
				ImGui::Combo("Combo", &n, "Yes\0No\0Maybe\0\0");
				ImGui::Checkbox("Check", &b);
				ImGui::EndMenu();
			}
			if (ImGui::BeginMenu("Colors"))
			{
				float sz = ImGui::GetTextLineHeight();
				for (int i = 0; i < ImGuiCol_COUNT; i++)
				{
					const char* name = ImGui::GetStyleColorName((ImGuiCol)i);
					ImVec2 p = ImGui::GetCursorScreenPos();
					ImGui::GetWindowDrawList()->AddRectFilled(p, ImVec2(p.x + sz, p.y + sz), ImGui::GetColorU32((ImGuiCol)i));
					ImGui::Dummy(ImVec2(sz, sz));
					ImGui::SameLine();
					ImGui::MenuItem(name);
				}
				ImGui::EndMenu();
			}
			if (ImGui::BeginMenu("Disabled", false)) // Disabled
			{
				IM_ASSERT(0);
			}
			
			if (ImGui::MenuItem("Quit", "Alt+F4")) { exit(0);  }
			ImGui::EndMenu();
		}
		if (ImGui::BeginMenu("View"))
		{
			if (ImGui::MenuItem("File menu", NULL, show_file_menu, true)) { show_file_menu = !show_file_menu; }
			if (ImGui::MenuItem("Legend", NULL, show_legend, true)) { show_legend = !show_legend; }
			ImGui::Separator();
			if (ImGui::MenuItem("Cut", "CTRL+X")) {}
			if (ImGui::MenuItem("Copy", "CTRL+C")) {}
			if (ImGui::MenuItem("Paste", "CTRL+V")) {}
			ImGui::EndMenu();
		}
		ImGui::EndMainMenuBar();
	}
}


// -----------------------------------------------------------------------------
static float padding_general = 5;
static float padding_top = 38;
static float padding_left = 5;
static float padding_right = 5;
static float menu_y = 190;
static float height_file_menu = 200;
static float height_points_menu = 200;
static float height_mesh_menu = 200;
static float height_histogram = 200;
static float height_legend = 310;
static float clicking_menu_height = 450;
static float viewer_menu_height = 350;
static float menu_width = 300;

static const ImGuiWindowFlags main_window_flags = ImGuiWindowFlags_NoScrollbar |
	ImGuiWindowFlags_NoResize |
	ImGuiWindowFlags_NoCollapse |
	ImGuiWindowFlags_NoNav;

void UIState::draw_file_menu(int x, int y) {
	ImGui::SetNextWindowPos(ImVec2(x,y), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, height_file_menu), ImGuiSetCond_FirstUseEver);

	ImGui::Begin("Cellogram", &show_file_menu, main_window_flags );

	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;


	if (ImGui::Button("Load Image", ImVec2((w - p), 0))) {
		std::string fname = FileDialog::openFileName(DATA_DIR, { "*.png", "*.tif", "*.tiff" });
		if (!fname.empty()) {
			load_image(fname);

			show_image = true;

			// update UI
			viewer_control();
		}
	}

	if (ImGui::Button("Save##Points", ImVec2((w - p) / 2.f, 0))) {
		save();
	}

	ImGui::SameLine();

	if (ImGui::Button("Load##Points", ImVec2((w - p) / 2.f, 0))) {
		load();
	}

	ImGui::End();
}

void UIState::draw_points_menu(int x, int y) {
	ImGui::SetNextWindowPos(ImVec2(x, y), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, height_points_menu), ImGuiSetCond_FirstUseEver);

	ImGui::Begin("Points", &show_file_menu, main_window_flags);

	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;
	//-------- Points ---------

	ImGui::InputFloat("Sigma", &state.sigma, 0.1, 0, 2);

	if (ImGui::Button("Detection", ImVec2((w - p), 0))) {
		detect_vertices();
	}

	ImGui::End();
}

void UIState::draw_mesh_menu(int x, int y)
{
	ImGui::SetNextWindowPos(ImVec2(x, y), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, 400), ImGuiSetCond_FirstUseEver);

	ImGui::Begin("Mesh", &show_file_menu, main_window_flags);

	//-------- Mesh ----------

	float w = ImGui::GetContentRegionAvailWidth();
	float p = ImGui::GetStyle().FramePadding.x;

	ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.40f);
	// Lloyds relaxation panel

	// disable if points are not detected
	if (state.mesh.points.size() == 0) push_disabled();

	ImGui::InputInt("Num Iter", &state.lloyd_iterations);

	if (ImGui::SliderFloat("t", &t, 0, 1)) {
		viewer_control();
	}

	ImGui::PopItemWidth();

#ifdef WITH_UNTANGLER
	if (ImGui::Button("Untangle", ImVec2((w - p), 0))) {
		state.untangle();
		t = 1;
		mesh_color.resize(0, 0);
		viewer_control();
	}

	ImGui::Separator();
#endif

	if (ImGui::Button("Lloyd", ImVec2((w - p), 0))) {
		state.relax_with_lloyd();
		t = 1;
		mesh_color.resize(0, 0);
		viewer_control();
	}

	if (ImGui::SliderFloat("energy", &state.energy_variation_from_mean, 0, 5)) {

		state.detect_bad_regions();
		state.check_regions();

		create_region_label();

		viewer_control();

	}

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

	if (state.mesh.points.size() == 0) pop_disabled();

	// disable if regions are not availabe
	if (state.regions.size() == 0) push_disabled();

	if (ImGui::Button("solve regions", ImVec2((w - p), 0))) {
		state.resolve_regions();

		create_region_label();

		viewer_control();
	}

	if (ImGui::Button("Ultimate relaxation", ImVec2((w - p), 0))) {
		state.final_relax();
		t = 1;
		create_region_label();

		viewer_control();
	}

	if (state.regions.size() == 0) pop_disabled();

	ImGui::End();
}

void UIState::draw_histogram(int x, int y)
{
	ImGui::SetNextWindowPos(ImVec2(x, y), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, height_histogram), ImGuiSetCond_FirstUseEver);
	ImGui::Begin("Histogram", &show_histogram, main_window_flags);
	const float hist_w = ImGui::GetWindowWidth() * 0.80f -2;
	ImGui::PushItemWidth(hist_w+2);
	// Image Scaling
	const int n_bins = 50;
	Eigen::Matrix<float, n_bins + 1, 1> hist;
	hist.setZero();

	for (long i = 0; i < state.img.size(); ++i)
	{
		hist(state.img(i)*n_bins)++;
	}


	static float min_img = 0;
	static float max_img = 1;

	auto pos = ImGui::GetWindowPos();
	int startX = pos.x + 10;
	int startY = pos.y + 42;

	ImGui::PushItemFlag(ImGuiItemFlags_Disabled, true);
	ImGui::PlotHistogram("", hist.data(), hist.size(), 0, NULL, 0.0f, hist.maxCoeff(), ImVec2(0, 80));

	ImDrawList* draw_list = ImGui::GetWindowDrawList();
	draw_list->PushClipRectFullScreen();
	draw_list->AddLine(ImVec2(startX + hist_w * min_img, startY), ImVec2(startX + hist_w * min_img, startY + 75), IM_COL32(0, 255, 0, 255), 2.0f);
	draw_list->AddLine(ImVec2(startX + hist_w * max_img, startY), ImVec2(startX + hist_w * max_img, startY + 75), IM_COL32(0, 255, 0, 255), 2.0f);
	draw_list->PopClipRect();

	ImGui::PopItemFlag();

	static const auto clamping = [](const double x) {
		if (x > 0) {
			if (x < 1)
				return x;
			else
				return 1.;
		}
		else
			return 0.0;
	};

	if (ImGui::SliderFloat("min", &min_img, 0.f, 1.f))
	{
		Eigen::MatrixXd tmp = (state.img.array() - min_img) / (max_img - min_img);
		tmp = tmp.unaryExpr(clamping);

		texture = (tmp.array() * 255).cast<unsigned char>();

		viewer_control();
	}
	if (ImGui::SliderFloat("max", &max_img, 0.f, 1.f))
	{
		Eigen::MatrixXd tmp = (state.img.array() - min_img) / (max_img - min_img);
		tmp = tmp.unaryExpr(clamping);
		texture = (tmp.array() * 255).cast<unsigned char>();

		viewer_control();
	}

	ImGui::PopItemWidth();
	ImGui::End();
}

void UIState::draw_legend(int x, int y) {
	ImGui::SetNextWindowPos(ImVec2(x, y), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, height_legend), ImGuiSetCond_FirstUseEver);
	if (!ImGui::Begin("Legend", &show_legend, main_window_flags))
	{
		ImGui::End();
		return;
	}

	draw_legend_item(46, 204, 113, "Ok");
	draw_legend_item(155, 89, 182, "Too Many Vertices");
	draw_legend_item(241, 196, 15, "Too Few Vertices");
	draw_legend_item(41, 128, 185, "Region Too Large");
	draw_legend_item(192, 57, 43, "No Solution");
	draw_legend_item(149, 165, 166, "Not Properly Closed");
	draw_legend_item(52, 73, 94, "Other");

	ImGui::End();
}

void UIState::draw_custom_window() {
	build_menu_bar();

	int x = padding_left, y = padding_top;

	// Menu on left
	if (show_file_menu)
	{
		draw_file_menu(x, y);
		y += height_file_menu + padding_general;
	}
	if (show_points_menu)
	{
		draw_points_menu(x, y);
		y += height_points_menu + padding_general;
	}
	if (show_mesh_menu)
	{
		draw_mesh_menu(x, y);
		y += height_mesh_menu + padding_general;
	}

	// Menu on right
	auto canvas = ImGui::GetIO().DisplaySize;
	x = canvas.x - padding_left - menu_width;
	y = padding_top;
	if (show_histogram)
	{
		draw_histogram(x, y);
		y += height_histogram + padding_general;
	}
	if (show_legend)
	{
		draw_legend(x, y);
		y += height_legend + padding_general;
	}
	//------------------------------------------------//
	//------------- Viewer Options--------------------//
	//------------------------------------------------//
	/*ImGui::SetNextWindowPos(ImVec2(menu_offset, main_menu_height + 2* menu_offset), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, viewer_menu_height), ImGuiSetCond_FirstUseEver);
	ImGui::Begin(
		"View settings", nullptr,
		ImGuiWindowFlags_NoSavedSettings
	);

	if (ImGui::ColorEdit4("Mesh color", points_data().line_color.data(),
		ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel)) {
		viewer_control();
	}

	if (ImGui::Checkbox("Show hull", &show_hull)) {
		viewer_control();
	}
	if (ImGui::Checkbox("Show points", &show_points)) {
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
	}*/


	//------------------------------------------------//
	//------------- Clicking Menu --------------------//
	//------------------------------------------------//
	//ImGui::SetNextWindowPos(ImVec2(menu_offset, main_menu_height + viewer_menu_height + 3 * menu_offset), ImGuiSetCond_FirstUseEver);
	//ImGui::SetNextWindowSize(ImVec2(menu_width, clicking_menu_height), ImGuiSetCond_FirstUseEver);
	//ImGui::Begin(
	//	"Clicking", nullptr,
	//	ImGuiWindowFlags_NoSavedSettings
	//);
	//if (state.regions.size() == 0) push_disabled();
	//if (ImGui::Button("Select Region", ImVec2(-1, 0))) {
	//	create_region_label();
	//	select_region = true;
	//	show_selected_region = true;
	//}
	//if (state.regions.size() == 0) pop_disabled();

	//int n_was_selected = selected_region;
	//if (n_was_selected < 0) push_disabled();
	//if (ImGui::Button("Grow Selected", ImVec2(-1, 0))) {
	//	state.grow_region(selected_region);
	//	create_region_label();
	//	viewer_control();
	//}
	//if (ImGui::Button("Solve Selected", ImVec2(-1, 0))) {
	//	state.resolve_region(selected_region);

	//	selected_region = -1;
	//	current_region_status = "";
	//	create_region_label();
	//	viewer_control();
	//}
	//if (n_was_selected < 0) pop_disabled();

	//bool was_delete = delete_vertex;
	//if (was_delete) push_selected();
	//if (ImGui::Button("Delete Vertex", ImVec2(-1, 0))) {
	//	add_vertex = false;
	//	delete_vertex = !delete_vertex;
	//}
	//if (was_delete) pop_selected();

	//bool was_add = add_vertex;
	//if (was_add) push_selected();
	//if (ImGui::Button("Add Vertex", ImVec2(-1, 0))) {
	//	delete_vertex = false;
	//	add_vertex = !add_vertex;
	//}
	//if (was_add) pop_selected();


	//if (ImGui::Button("Mark good", ImVec2(-1, 0))) {
	//	// select vertices and mark them as good permanently
	//	make_vertex_good = true;
	//	viewer_control();
	//}
	//if (ImGui::Button("Mark bad", ImVec2(-1, 0))) {
	//	// select vertices and mark them as good permanently
	//	make_vertex_bad = true;
	//	viewer_control();
	//}

	//bool was_moved = move_vertex;
	//if (was_moved) push_selected();
	//if (ImGui::Button("Move vertex", ImVec2(-1, 0))) {
	//	// drag a single vertex to a new starting position
	//	move_vertex = !move_vertex;
	//	viewer_control();
	//}
	//if (was_moved) pop_selected();


	//if (ImGui::Button("Color code", ImVec2(-1, 0))) {
	//	// select vertices and mark them as good permanently
	//	color_code = !color_code;
	//	viewer_control();
	//}

	//if (selected_region >= 0) {
	//	ImGui::PushItemWidth(ImGui::GetWindowWidth());
	//	int nVi = state.regions[selected_region].region_interior.size();
	//	int nVtotal = state.regions[selected_region].region_boundary.size() + nVi;
	//	int nTri = state.regions[selected_region].region_triangles.size();
	//	ImGui::LabelText("", "Region %d ", selected_region);
	//	if (state.regions[selected_region].points_delta != 0)
	//		ImGui::LabelText("", "%s (%i) ", current_region_status.c_str(), state.regions[selected_region].points_delta);
	//	else
	//		ImGui::LabelText("", "%s ", current_region_status.c_str());
	//	ImGui::LabelText("", "#V: %i (%i) ", nVtotal, nVi);
	//	ImGui::LabelText("", "#F: %i ", nTri);
	//	ImGui::PopItemWidth();
	//}



	if (delete_vertex || add_vertex)
	{
		//Cross hair
		ImGui::SetNextWindowPos(ImVec2(-100, -100), ImGuiSetCond_Always);
		ImGui::Begin("mouse_layer");
		ImVec2 p = ImGui::GetIO().MousePos;
		ImDrawList* draw_list = ImGui::GetWindowDrawList();
		draw_list->PushClipRectFullScreen();
		draw_list->AddLine(ImVec2(p.x - 50, p.y), ImVec2(p.x + 50, p.y), IM_COL32(delete_vertex ? 255 : 0, add_vertex ? 255 : 0, 0, 255), 2.0f);
		draw_list->AddLine(ImVec2(p.x, p.y - 50), ImVec2(p.x, p.y + 50), IM_COL32(delete_vertex ? 255 : 0, add_vertex ? 255 : 0, 0, 255), 2.0f);
		draw_list->PopClipRect();

		ImGui::GetIO().MouseDrawCursor = true;
		ImGui::SetMouseCursor(-1);
		ImGui::End();
	}
	else
	{
		ImGui::GetIO().MouseDrawCursor = false;
		ImGui::SetMouseCursor(0);
	}


}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

} // namespace cellogram
