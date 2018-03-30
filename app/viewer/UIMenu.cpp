////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include "FileDialog.h"
#include <cellogram/voronoi.h>
#include <cellogram/laplace_energy.h>
#include <cellogram/tri2hex.h>
#include <cellogram/vertex.h>
#include <cellogram/region_grow.h>
#include <cellogram/vertex_degree.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
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
	ImGui::SetNextWindowPos(ImVec2(190.f * menu_scaling(), 0), ImGuiSetCond_FirstUseEver);
	ImGui::SetNextWindowSize(ImVec2(300, 600), ImGuiSetCond_FirstUseEver);
	ImGui::Begin(
		"Debug", nullptr,
		ImGuiWindowFlags_NoSavedSettings
	);
	// Lloyds relaxation panel
	if (ImGui::CollapsingHeader("Lloyds Relaxation", ImGuiTreeNodeFlags_DefaultOpen)) {
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		if (ImGui::Button("Lloyd", ImVec2((w - p) / 2.f, 0))) {
			lloyd_relaxation(state.points, state.boundary, 1, state.hull_vertices, state.hull_faces);
			points_data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
			compute_triangulation();
		}
		//ImGui::SameLine(0, p);
		ImGui::Checkbox("Animate", &continuous_lloyd);
		static float t = 0;
		if (ImGui::SliderFloat("t", &t, 0, 1)) {
			Eigen::MatrixXd V = t * state.points + (1 - t) * state.detected;
			points_data().set_mesh(V, state.triangles);
		}
	}

	// Laplace energy panel
	if (ImGui::CollapsingHeader("Laplace Energy", ImGuiTreeNodeFlags_DefaultOpen)) {
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		if (ImGui::Button("Current pos", ImVec2((w - p) / 2.f, 0))) {
			laplace_energy(state.points, state.triangles, state.current_laplace_energy);
			Eigen::MatrixXd C;
			igl::parula(state.current_laplace_energy, true, C);

			// Plot the mesh with pseudocolors
			points_data().clear();
			points_data().set_mesh(state.points, state.triangles);
			points_data().set_colors(C);
			fix_color(points_data());

		}
		ImGui::SameLine(0, p);
		if (ImGui::Button("Original pos", ImVec2((w - p) / 2.f, 0))) {
			laplace_energy(state.detected, state.triangles, state.original_laplace_energy);
			Eigen::MatrixXd C;
			igl::parula(state.original_laplace_energy, true, C);

			// Plot the mesh with pseudocolors
			points_data().clear();
			points_data().set_mesh(state.points, state.triangles);
			points_data().set_colors(C);
			fix_color(points_data());
		}
	}
	// Image Control Panel
	if (ImGui::CollapsingHeader("Image", ImGuiTreeNodeFlags_DefaultOpen)) {
		if (ImGui::Button("Load Image")) {
			std::string fname = FileDialog::openFileName(DATA_DIR, { "*.png" });
			if (!fname.empty()) {
				load_image(fname);
				img_data().show_texture = true;
				// now that image is loaded, menu for image handling can be called
				image_loaded = true;
				fix_color(img_data());
			}
		}

		// Generate new submenu to handle the image
		/*
		It should be able to
		- turn on and off the image
		- change the contrast
		-
		*/
		if (image_loaded) {
			ImGui::Checkbox("Show image", &(img_data().show_faces));
			// Brightness
			ImGui::PushItemWidth(80 * menu_scaling());
			ImGui::DragFloat("Brightness", &(viewer.core.camera_zoom), 0.05f, 0.1f, 20.0f);
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
	if (ImGui::CollapsingHeader("Vertices", ImGuiTreeNodeFlags_DefaultOpen)) {
		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;
		if (ImGui::Button("Add node", ImVec2((w - p) / 2.f, 0))) {
			add_vertex(state.points);
		}
		ImGui::SameLine(0, p);
		if (ImGui::Button("Delete Node", ImVec2((w - p) / 2.f, 0))) {
			delete_vertex(state.points);
		}
		ImGui::ColorEdit4("Node color", points_data().line_color.data(),
			ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
	}
	// Mesh
	/* This menu will include:
	- show face
	- changing color of mesh
	*/
	if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen)) {
		ImGui::Checkbox("Mesh Fill", &(points_data().show_faces));
		ImGui::ColorEdit4("Mesh color", points_data().line_color.data(),
			ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
	}
	
	// Button to call any function for testing
	if (ImGui::CollapsingHeader("Test foo", ImGuiTreeNodeFlags_DefaultOpen)) {
		if (ImGui::Button("Test function")) {

			// Calculate the graph adjancency
			state.adjacency_list.clear();
			tri2hex(state.triangles, state.adjacency_list);

			// Calculate the laplacian energy with respect to the original positions
			laplace_energy(state.points, state.triangles, state.current_laplace_energy);

			// Find the degree of each vertex
			Eigen::VectorXi degree;
			vertex_degree(state.triangles, degree );

			// Determine whether each vertex passes the criterium for bad mesh
			double avg = state.current_laplace_energy.mean();
			std::vector<bool> crit_pass(state.points.rows(), false);
			for (int i = 0; i < state.points.rows(); i++)
			{
				if (state.current_laplace_energy(i) > 0.8*avg || degree(i) != 6)
				//if (degree(i) != 6)
				{
					crit_pass[i] = true;
				}
			}

			// Find connected regions where the criterium was not passed
			Eigen::VectorXi region;
			region_grow(state.adjacency_list, crit_pass, region);

			//std::cout << region.transpose() << std::endl;

			// Find edges of connected region
			std::vector<std::vector<int>> region_edges;
			region_bounding(state.points, state.triangles, region, region_edges);
			/*
			std::cout << "Points\n" << state.points.transpose() << "\nBoundary\n" << std::endl;
			for (int j = 0; j < region_edges[0].size(); j++) {
				std::cout << region_edges[0][j] << " ";
			}*/
			

			// Color mesh according to 
			Eigen::MatrixXd C;
			Eigen::VectorXd regionD = region.cast<double>();
			igl::jet(regionD, true, C);
			points_data().clear();
			points_data().set_mesh(state.points, state.triangles);
			points_data().set_colors(C);
			fix_color(points_data());

			// Line boundary
			Eigen::MatrixXd bad_P1(state.points.rows(), 3);
			Eigen::MatrixXd bad_P2(state.points.rows(), 3);
			int k = 0;
			for (int i = 0; i < region_edges.size(); i++)
			{
				int n = region_edges[i].size();
				for (int j = 0; j < region_edges[i].size(); j++) //
				{
					bad_P1.row(k) = state.points.row(region_edges[i][j]);
					bad_P2.row(k) = state.points.row(region_edges[i][(j + 1) % n]);
					k++;
				}
			}
			
			bad_P1.conservativeResize(k, 3);
			bad_P2.conservativeResize(k, 3);

			bad_region_data().clear();
			bad_region_data().add_edges(bad_P1, bad_P2, Eigen::RowVector3d(0, 0, 0));
			
			//std::cout << "Boundary\n" << bad_P1.transpose() << std::endl;
			//std::cout << "Points\n" << state.points.transpose() << std::endl;

			// Draw filled polygon
			Eigen::MatrixXd bad_region_vertices;
			Eigen::MatrixXi bad_region_faces;
			//std::cout << bad_P1.transpose() << std::endl;
			triangulate_polygon(bad_P1, bad_region_vertices, bad_region_faces);
			bad_region_data().set_mesh(bad_region_vertices, bad_region_faces);
			
			// Set viewer options
			bad_region_data().set_colors(Eigen::RowVector3d(52, 152, 219) / 255.0);
			bad_region_data().show_faces = false;
			bad_region_data().show_lines = false;
			bad_region_data().shininess = 0;
			bad_region_data().line_width = 4.0;
			
		}
	}

	ImGui::End();
}

// -----------------------------------------------------------------------------

bool UIState::pre_draw() {
	ImGuiMenu::pre_draw();

	// perform lloyd if necessary
	if (continuous_lloyd) {
		lloyd_relaxation(state.points, state.boundary, 1, state.hull_vertices, state.hull_faces);
		points_data().set_points(state.points, Eigen::RowVector3d(1, 0, 0));
		compute_triangulation();
	}

	return false;
}

// -----------------------------------------------------------------------------

} // namespace cellogram
