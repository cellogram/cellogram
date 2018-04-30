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

// -----------------------------------------------------------------------------

void UIState::draw_viewer_menu() {
	// Workspace
	// if (ImGui::CollapsingHeader("Workspace", ImGuiTreeNodeFlags_DefaultOpen)) {
	// 	float w = ImGui::GetContentRegionAvailWidth();
	// 	float p = ImGui::GetStyle().FramePadding.x;
	// 	if (ImGui::Button("Load##Workspace", ImVec2((w-p)/2.f, 0))) {
	// 		viewer.load_scene();
	// 	}
	// 	ImGui::SameLine(0, p);
	// 	if (ImGui::Button("Save##Workspace", ImVec2((w-p)/2.f, 0))) {
	// 		viewer.save_scene();
	// 	}
	// }

	// Mesh
	if (ImGui::CollapsingHeader("Points", ImGuiTreeNodeFlags_DefaultOpen)) {
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		/*if (ImGui::Button("Load##Points", ImVec2((w-p)/2.f, 0))) {
			std::string fname = FileDialog::openFileName(DATA_DIR, {"*.xyz"});
			if (!fname.empty()) { load(fname); }
		}*/
		//ImGui::SameLine(0, p);
		if (ImGui::Button("Save##Points", ImVec2((w-p)/2.f, 0))) {
			std::string fname = FileDialog::saveFileName(DATA_DIR, {"*.xyz"});
			if (!fname.empty()) { save(fname); }
		}
	}

	// Viewing options
	if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen)) {
		if (ImGui::Button("Center object", ImVec2(-1, 0))) {
			viewer.core.align_camera_center(viewer.data().V, viewer.data().F);
		}
		if (ImGui::Button("Snap canonical view", ImVec2(-1, 0))) {
			viewer.snap_to_canonical_quaternion();
		}

		// Zoom
		ImGui::PushItemWidth(80 * menu_scaling());
		ImGui::DragFloat("Zoom", &(viewer.core.camera_zoom), 0.05f, 0.1f, 20.0f);

		// Select rotation type
		static int rotation_type = static_cast<int>(viewer.core.rotation_type);
		static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
		static bool orthographic = true;
		if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0")) {
			using RT = igl::opengl::ViewerCore::RotationType;
			auto new_type = static_cast<RT>(rotation_type);
			if (new_type != viewer.core.rotation_type) {
				if (new_type == RT::ROTATION_TYPE_NO_ROTATION) {
					trackball_angle = viewer.core.trackball_angle;
					orthographic = viewer.core.orthographic;
					viewer.core.trackball_angle = Eigen::Quaternionf::Identity();
					viewer.core.orthographic = true;
				} else if (viewer.core.rotation_type == RT::ROTATION_TYPE_NO_ROTATION) {
					viewer.core.trackball_angle = trackball_angle;
					viewer.core.orthographic = orthographic;
				}
				viewer.core.set_rotation_type(new_type);
			}
		}

		// Orthographic view
		ImGui::Checkbox("Orthographic view", &(viewer.core.orthographic));
		ImGui::PopItemWidth();
	}

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

	// Overlays
	if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen)) {
		ImGui::Checkbox("Wireframe", &(viewer.data().show_lines));
		ImGui::Checkbox("Fill", &(viewer.data().show_faces));
		ImGui::Checkbox("Show vertex labels", &(viewer.data().show_vertid));
		ImGui::Checkbox("Show faces labels", &(viewer.data().show_faceid));
	}
}

// -----------------------------------------------------------------------------
static float menu_y = 190;
static float menu_height = 750;
static float menu_width = 300;
void UIState::draw_custom_window() {
	ImGui::SetNextWindowPos(ImVec2(menu_y * menu_scaling(), 0), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, menu_height), ImGuiSetCond_FirstUseEver);
	ImGui::Begin(
		"Algorithm", nullptr,
		ImGuiWindowFlags_NoSavedSettings
	);
	// Lloyds relaxation panel
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;

		ImGui::InputInt("Num Iter", &state.lloyd_iterations);
		/*if (ImGui::Button("Lloyd", ImVec2((w - p) / 2.f, 0))) {
			state.relax_with_lloyd();
			t = 1;
			mesh_color.resize(0, 0);
			viewer_control();
		}*/
		
	//}

	if (ImGui::SliderFloat("t", &t, 0, 1)) {
		viewer_control();
	}

	// Laplace energy panel
	if (ImGui::CollapsingHeader("Laplace Energy", ImGuiTreeNodeFlags_DefaultOpen)) {
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		if (ImGui::Button("Current pos", ImVec2((w - p) / 2.f, 0))) {
			Eigen::VectorXd current_laplace_energy;
			laplace_energy(state.mesh.points, state.mesh.triangles, current_laplace_energy);
			igl::parula(current_laplace_energy, true, mesh_color);

			// update UI
			viewer_control();

		}
		ImGui::SameLine(0, p);
		if (ImGui::Button("Original pos", ImVec2((w - p) / 2.f, 0))) {
			Eigen::VectorXd original_laplace_energy;
			laplace_energy(state.mesh.detected, state.mesh.triangles, original_laplace_energy);
			igl::parula(original_laplace_energy, true, mesh_color);

			// update UI
			viewer_control();
		}
	}

	// Node control panel
	/* This menu will include:
		- adding/removing nodes
		- moving points
		- changing color of nodes
			. with meaning, e.g. valence
			. color picked by user
	*/
	/*if (ImGui::CollapsingHeader("Vertices", ImGuiTreeNodeFlags_DefaultOpen)) {
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		if (ImGui::Button("Add node", ImVec2((w - p) / 2.f, 0))) {

			//add_vertex();
		}
		ImGui::SameLine(0, p);
		if (ImGui::Button("Delete Node", ImVec2((w - p) / 2.f, 0))) {
			//delete_vertex(state.points);
		}
		if (ImGui::ColorEdit4("Vertex color", vertex_color.data(),
			ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel))
		{
			viewer_control();
		}
	}*/
	// Mesh
	/* This menu will include:
	- show face
	- changing color of mesh
	*/
	if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen)) {
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		if (ImGui::ColorEdit4("Mesh color", points_data().line_color.data(),
			ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel)) {
			viewer_control();
		}
	}
	
	if (ImGui::CollapsingHeader("Phases", ImGuiTreeNodeFlags_DefaultOpen)) {
		// Lloyds relaxation panel
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;

		if (ImGui::Button("Untangle", ImVec2((w - p) , 0))) {
			state.untangle();
			t = 1;
			mesh_color.resize(0, 0);
			viewer_control();
		}

		ImGui::Separator();

		if (ImGui::Button("Load Image")) {
			std::string fname = FileDialog::openFileName(DATA_DIR, { "*.png" });
			if (!fname.empty()) {
				load_image(fname);

				show_image = true;

				// update UI
				viewer_control();
			}
		}

		if (ImGui::Button("Detection")) {
			
			detect_vertices();
		}

		if (ImGui::Button("Lloyd", ImVec2((w - p) , 0))) {
			state.relax_with_lloyd();
			t = 1;
			mesh_color.resize(0, 0);
			viewer_control();
		}


		if (ImGui::Button("build regions", ImVec2((w - p), 0))) {
			state.detect_bad_regions();
			
			show_points = false;
			show_bad_regions = true;

			igl::jet(create_region_label(), true, mesh_color);
			
			viewer_control();
		}

		if (ImGui::Button("check regions", ImVec2((w - p) , 0))) {
			state.check_regions();
		}

		if (ImGui::Button("fix regions", ImVec2((w - p), 0))) {
			state.fix_regions();

			igl::jet(create_region_label(), true, mesh_color);

			viewer_control();
		}

		if (ImGui::Button("solve regions", ImVec2((w - p), 0))) {
			state.resolve_regions();			

			igl::jet(create_region_label(), true, mesh_color);

			viewer_control();
		}

		if (ImGui::Button("grow regions", ImVec2((w - p) , 0))) {
			state.grow_regions();			

			igl::jet(create_region_label(), true, mesh_color);

			viewer_control();
		}

		if (ImGui::Button("Ultimate relaxation", ImVec2((w - p), 0))) {
			state.final_relax();

			igl::jet(create_region_label(), true, mesh_color);

			viewer_control();
		}
	}

	ImGui::End();


	// Clicking Menu
	ImGui::SetNextWindowPos(ImVec2(0 * menu_scaling(), menu_height + 5), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, 350), ImGuiSetCond_FirstUseEver);
	ImGui::Begin(
		"Clicking", nullptr,
		ImGuiWindowFlags_NoSavedSettings
	);

	if (ImGui::Button("Select Region", ImVec2(-1, 0))) {
		igl::jet(create_region_label(), true, mesh_color);
		select_region = true;
	}
	if (ImGui::Button("Grow Selected", ImVec2(-1, 0))) {
		if (selected_region < 0) return;
		state.grow_region(selected_region);
		igl::jet(create_region_label(), true, mesh_color);
		viewer_control();
	}
	if (ImGui::Button("Solve Selected", ImVec2(-1, 0))) {
		state.resolve_region(selected_region);

		selected_region = -1;
		current_region_status = "";
		igl::jet(create_region_label(), true, mesh_color);
		viewer_control();
	}
	if (ImGui::Button("Delete Vertex", ImVec2(-1, 0))) {
		delete_vertex = true;
		
	}
	if (ImGui::Button("Add Vertex", ImVec2(-1, 0))) {
		add_vertex = true;
	}
	ImGui::InputInt("Param", &selected_param);
	if (ImGui::Button("Color Code", ImVec2(-1, 0))) {
		// determine color map for interior vertices
		color_code = true;
		show_points = true;
		viewer_control();
	}
	if (ImGui::Button("Split region", ImVec2(-1, 0))) {
		// select vertices and mark them as good permenantly
		make_vertex_good = true;
		viewer_control();
	}

	if (selected_region >= 0) {
		ImGui::LabelText("", "Region %d", selected_region);
		ImGui::LabelText("", "%s", current_region_status.c_str());
	}


	ImGui::End();



	// Displaying Options
	ImGui::SetNextWindowPos(ImVec2(menu_y * menu_scaling(), menu_height + 5), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(menu_width, 300), ImGuiSetCond_FirstUseEver);
	ImGui::Begin(
		"View settings", nullptr,
		ImGuiWindowFlags_NoSavedSettings
	);

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
	ImGui::End();
}

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------

} // namespace cellogram
