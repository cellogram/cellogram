#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "cellogram/State.h"
#include <Eigen/Dense>
#include <igl/colormap.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>

#include <imgui/imgui.h>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

typedef struct UIsize_t {
	bool resize;
    int windowWidth;
    int windowHeight;
    int leftWidth;
	int rightWidth;
    int bottomHeight;
    int imageViewerHeight;
    double mainMenuHeight;

	UIsize_t() : resize(true), windowWidth(1600), windowHeight(900), 
				 leftWidth(300), rightWidth(300), bottomHeight(150), 
				 imageViewerHeight(350) {}
} UTsize_t;


typedef struct ImageViewer_t {
	int imageViewerType;
	int sliceToShow;
	float visual_z_mult;  // y-axis visualization multiplier (z-stack may be too small)
	float deformScale;  // mesh deformation percent
	float darkenScale;

	ImageViewer_t() : imageViewerType(0), sliceToShow(0), visual_z_mult(5.0), deformScale(0.0), darkenScale(0.5) {}
} ImageViewer_t;


typedef struct ColorUI_t {
    Eigen::RowVector3f mesh_vertex_color;
    Eigen::MatrixXd mesh_fill_color;
    Eigen::RowVector3f mesh_line_color;
    Eigen::RowVector3d match_line_color;

    Eigen::MatrixXd warn_snd_derivative;
    Eigen::MatrixXd warn_invalid_energy;

    Eigen::MatrixXd radius_point_color;
    Eigen::MatrixXd ref_point_color;  // reference axis color
    Eigen::Vector4f background_color;

    ColorUI_t() {
        mesh_vertex_color = Eigen::RowVector3f(0, 152, 255)/255.0;
        mesh_fill_color.resize(1, 3);
        mesh_fill_color << 0.8, 0.8, 0.85;
        mesh_line_color = Eigen::RowVector3f(3, 21, 20)/255.0;
        match_line_color = Eigen::RowVector3d(194./255., 83./255., 57./255.);

        warn_snd_derivative.resize(1, 3);
        warn_snd_derivative << 250./255., 195./255., 92./255.;
        warn_invalid_energy.resize(1, 3);
        warn_invalid_energy << 250./255., 55./255., 41./255.;

        radius_point_color.resize(1, 3);
        radius_point_color << 255./255., 133./255., 44./255.;
        ref_point_color.resize(1, 3);
        ref_point_color << 0., 0., 0.3;
        background_color = Eigen::Vector4f(0.7, 0.7, 0.75, 1.0);
    }
} ColorUI_t;

// -----------------------------------------------------------------------------

class UIState : public igl::opengl::glfw::imgui::ImGuiMenu {
private:
	UIState();

public:
	typedef igl::opengl::glfw::imgui::ImGuiMenu super;
	static UIState &ui_state();

	virtual ~UIState() = default;

public:
	// Viewer
	igl::opengl::glfw::Viewer viewer;

	// Scene
	State &state;

	// logger
	std::ostringstream oss;

	// Multiple meshes ids
	int hull_id;
	int points_id;
	int image_id;
	int bad_region_id;
	int matching_id;
	int selected_id;
	int physical_id;
	int visual_id;
    int warn_id;

	int selected_region = -1;
	int selected_param = 0;
	int dragging_id = -1;
	// UI options
	// double foo;

	Eigen::MatrixXi img_F;
	Eigen::MatrixXd img_V;
	Eigen::MatrixXf hist;

	// Display flags
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> texture;
    ColorUI_t colorUI;

	bool show_mesh = true;
	bool show_hull = false;
	//bool image_loaded = false;
	bool show_points = false;
	bool show_mesh_fill = true;
	bool show_image = false;
	bool show_matching = false;
	bool show_bad_regions = false;
	bool color_code = false;
	bool show_selected_region = true;
	bool analysis_mode = false;
	bool show_traction_forces = true;
	bool show_smoothed_results = false;

	// 3d visualizer
	enum class Mesh3DAttribute : int {
		NONE,
		X_DISP,
		Y_DISP,
		Z_DISP,
		NORM_DISP,
	};
	Mesh3DAttribute view_mode_3d = Mesh3DAttribute::NONE;

	// Clicking flags
	bool select_region = false;
	bool add_vertex = false;
	bool delete_vertex = false;
	bool make_vertex_good = false;
	bool make_vertex_bad = false;
	bool move_vertex = false;
	int split_region = -1;
	Eigen::Vector2i split_end_points;

	// Image
	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> img;

	std::string current_region_status;

	std::string save_dir = "";
	bool data_available = false;
public:
	virtual void init(igl::opengl::glfw::Viewer *_viewer) override;
	void initialize();

	void launch();

	igl::opengl::ViewerData & mesh_by_id(int id);

	bool load();
	//bool load_param(std::string name);

	bool save();
	bool save_as(const std::string & save_as_dir);

	void post_resize(int w, int h) override;
	virtual bool mouse_down(int button, int modifier) override;
	virtual bool mouse_move(int button, int modifier) override;
	virtual bool mouse_up(int button, int modifier) override;
	// virtual bool mouse_scroll(float delta_y) override;

	void load_image(std::string name);
	void detect_vertices();
	void display_image();
	void compute_hull();
	void clean_hull();
	void compute_triangulation();
	void compute_histogram();

	void export_region();
	void reset_viewer();
	void deselect_all_buttons();
public:
	igl::opengl::ViewerData & points_data() { return mesh_by_id(points_id); }
	igl::opengl::ViewerData & hull_data() { return mesh_by_id(hull_id); }
	igl::opengl::ViewerData & image_data() { return mesh_by_id(image_id); }
	igl::opengl::ViewerData & bad_region_data() { return mesh_by_id(bad_region_id); }
	igl::opengl::ViewerData & matching_data() { return mesh_by_id(matching_id); }
	igl::opengl::ViewerData & selected_data() { return mesh_by_id(selected_id); }
	igl::opengl::ViewerData & physical_data() { return mesh_by_id(physical_id); }
	igl::opengl::ViewerData & visual_data() { return mesh_by_id(visual_id); }
    igl::opengl::ViewerData & warn_data() { return mesh_by_id(warn_id); }

private:
	igl::ColorMapType cm = igl::ColorMapType::COLOR_MAP_TYPE_PARULA;
	double min_val = 0, max_val = 0;
	double real_min_val = 0, real_max_val = 0;
	bool override_ranges = false;
	std::string current_file_name = "";
	ImFont *icon_font;
	ImFont *icon_font_big;

	bool block_mouse_behavior(int button);
	void viewer_control();
	void viewer_control_2d();
	void viewer_control_3d();
	void fix_color(igl::opengl::ViewerData &data);
	void create_region_label();
	void build_region_edges(const Eigen::MatrixXd &pts, Eigen::MatrixXd &bad_P1, Eigen::MatrixXd &bad_P2, Eigen::MatrixXd &C);

	void phase_0();
	void phase_1();
	void phase_2();
	void phase_3();
	void phase_4();
	void phase_5();
	///////////////
	// UI Panels //
	///////////////

	// Menu bar
	float draw_menu_bar(); // returns menu bar height

	// Left panel
	void draw_left_panel(float ypos, float width);
	float draw_file_menu();
	void draw_points_menu();
	void draw_mesh_menu();
	void draw_depth_menu();
	void draw_analysis_menu();
	void draw_results_menu();

	// Right panel
	void draw_right_panel(float ypos, float width);
	void draw_histogram_menu();
	void draw_legend_menu();
	void draw_layer_menu();
	void draw_region_menu();

	// image viewer
	void draw_image_viewer_menu();
	// log window
	void draw_log_window();
    // editor
    void draw_editor_window();

	// Toggle windows
	bool show_left_panel = true;
	bool show_file_menu = true;
	bool show_points_menu = true;
	bool show_mesh_menu = true;
	bool show_analysis_menu = true;

	bool show_right_panel = true;
	bool show_histogram_menu = true;
	bool show_legend_menu = true;
	bool show_layer_menu = true;
	bool show_region_menu = true;

	bool show_imageViewer_menu = true;
	bool show_log_menu = true;
    bool show_editor_menu = false;

	UIsize_t UIsize;
	ImageViewer_t imgViewer;

	// visualization
	bool show_axisPoints = false;
    bool show_allIndex = false;
    bool show_radiusPoints = false;
	void DrawAxisDots();
    void DrawAllMarkerIdx();
    void DrawWarnViewer();
    void DrawRadiusDots();

public:
	// Menu stuff
	void draw_viewer_window() override;

	bool key_pressed(unsigned int unicode_key, int modifiers) override;
	bool key_up(int key, int modifiers) override;

	std::function<void(void)> next_call_back;
	bool has_next_callback = false;
	bool close_next = false;
	int run_next = 0;
	template<typename BtnFun, typename CallBackFun>
	void wait_window(const std::string & id, const std::string & msg, const char* icon, const BtnFun & button, const CallBackFun & call_back)
	{
		std::string title = "Please wait##" + id;
		if (button()) {
			ImGui::OpenPopup(title.c_str());
			run_next = 0;

			close_next = false;
			next_call_back = call_back;
			has_next_callback = true;
		}
		ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiCond_Always);
		ImGui::SetNextWindowSizeConstraints(ImVec2(200, 60), ImVec2(200, 60));
		if (ImGui::BeginPopupModal(title.c_str(), NULL, ImGuiWindowFlags_AlwaysAutoResize))
		{
			ImGui::Spacing();
			ImGui::PushFont(icon_font);
			ImGui::Text("%s", icon);
			ImGui::PopFont();

			ImGui::SameLine();

			ImGui::Text("%s", msg.c_str());

			if(close_next)
			{
				ImGui::CloseCurrentPopup();
				close_next = false;
			}
			ImGui::EndPopup();
		}
	}
};

// -----------------------------------------------------------------------------

} // namespace cellogram
