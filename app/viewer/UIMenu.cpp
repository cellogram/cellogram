////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include "FileDialog.h"
#include <cellogram/voronoi.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui_internal.h>
#include <imgui/imgui.h>
#include <algorithm>
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
		if (ImGui::Button("Load##Points", ImVec2((w-p)/2.f, 0))) {
			std::string fname = FileDialog::openFileName(DATA_DIR, {"*.xyz"});
			if (!fname.empty()) { load(fname); }
		}
		ImGui::SameLine(0, p);
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
	if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen)) {
		if (ImGui::Checkbox("Face-based", &(viewer.data().face_based))) {
			viewer.data().set_face_based(viewer.data().face_based);
		}
		ImGui::Checkbox("Show texture", &(viewer.data().show_texture));
		if (ImGui::Checkbox("Invert normals", &(viewer.data().invert_normals))) {
			viewer.data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
		}
		ImGui::Checkbox("Show overlay", &(viewer.data().show_overlay));
		ImGui::Checkbox("Show overlay depth", &(viewer.data().show_overlay_depth));
		ImGui::ColorEdit4("Background", viewer.core.background_color.data(),
				ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
		ImGui::ColorEdit4("Line color", viewer.data().line_color.data(),
				ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
		ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
		ImGui::DragFloat("Shininess", &(viewer.data().shininess), 0.05f, 0.0f, 100.0f);
		ImGui::PopItemWidth();
	}

	// Overlays
	if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen)) {
		ImGui::Checkbox("Wireframe", &(viewer.data().show_lines));
		ImGui::Checkbox("Fill", &(viewer.data().show_faces));
		ImGui::Checkbox("Show vertex labels", &(viewer.data().show_vertid));
		ImGui::Checkbox("Show faces labels", &(viewer.data().show_faceid));
	}
}

// -----------------------------------------------------------------------------

void UIState::draw_custom_window() {
	if (ImGui::Button("Lloyd")) {
		lloyd_relaxation(state.points, Eigen::VectorXi(1), 1, state.hull, state.hull_faces);
		points_data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
	}
}

// -----------------------------------------------------------------------------

} // namespace cellogram
