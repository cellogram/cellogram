////////////////////////////////////////////////////////////////////////////////
#include "cellogram/convex_hull.h"
#include "cellogram/voronoi.h"
#include "viewer/file_dialog.h"
#include <geogram/basic/file_system.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram_gfx/glup_viewer/glup_viewer.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
#include <geogram_gfx/mesh/mesh_gfx.h>
#include <algorithm>
////////////////////////////////////////////////////////////////////////////////

typedef GEO::vec2 vec2;

struct Scene {
	std::vector<bool> fixed;
	std::vector<vec2> detected;
	std::vector<vec2> reference;

	std::vector<vec2> border;
	GEO::Mesh border_mesh;
};

////////////////////////////////////////////////////////////////////////////////

namespace {

	// Load vertices from a mesh file (mesh can be a .obj, .ply, .xyz, etc.)
	std::vector<vec2> load_points(const std::string &filename) {
		GEO::Mesh M;
		if(!GEO::mesh_load(filename, M)) {
			return {};
		}
		std::vector<vec2> points;
		for (GEO::index_t i = 0; i < M.vertices.nb(); ++i) {
			GEO::vec3 p = M.vertices.point(i);
			points.push_back(vec2(p[0], p[1]));
		}
		return points;
	}

	void update_mesh(const std::vector<vec2> &points, GEO::Mesh &M) {
		for (GEO::index_t i = 0; i < M.vertices.nb(); ++i) {
			if (M.vertices.single_precision()) {
				std::copy_n(&points[i][0], 2, M.vertices.single_precision_point_ptr(i));
			} else {
				std::copy_n(&points[i][0], 2, M.vertices.point_ptr(i));
			}
		}
	}

} // anonymous namespace

////////////////////////////////////////////////////////////////////////////////

class CellogramApplication : public GEO::SimpleMeshApplication {
	Scene scene;

	bool show_detected_;
	bool show_reference_;
	bool show_border_;

public:
	CellogramApplication(int argc, char** argv, const std::string& usage)
		: SimpleMeshApplication(argc, argv, "")
		, show_detected_(true)
		, show_reference_(false)
		, show_border_(true)
	{
		name_ = "[Float] Cellogram Viewer";
	}

	virtual void init_graphics() override {
		SimpleMeshApplication::init_graphics();
		glup_viewer_disable(GLUP_VIEWER_BACKGROUND);
		glup_viewer_disable(GLUP_VIEWER_3D);
	}

	virtual void draw_load_menu() override {
		if(ImGui::MenuItem("Load detected...")) {
			std::string filename = FileDialog::openFileName(DATA_DIR);
			load(filename);
		}
	}

	virtual bool load(const std::string &filename) override {
		if(!GEO::FileSystem::is_file(filename)) {
			GEO::Logger::out("I/O") << "is not a file" << std::endl;
			return false;
		}
		scene.detected = load_points(filename);
		SimpleMeshApplication::load(filename);
		current_file_ = "";
		auto hull = cellogram::convex_hull(scene.detected);
		//std::reverse(hull.begin(), hull.end());
		scene.fixed.assign(scene.detected.size(), false);
		scene.border.clear();
		for (int x : hull) {
			scene.fixed[x] = true;
			scene.border.push_back(scene.detected[x]);
		}
		cellogram::triangulate_hull(scene.border, scene.border_mesh);
		return true;
	}

	virtual void draw_save_menu() override {
		if(ImGui::MenuItem("Save...")) {
			std::string filename = FileDialog::saveFileName();
		}
	}

	virtual void draw_viewer_properties() override {
		ImGui::Text("Display");
		ImGui::Checkbox("detected points", &show_detected_);
		ImGui::Checkbox("reference points", &show_reference_);
		ImGui::Checkbox("border", &show_border_);
		ImGui::Separator();
		ImGui::Text("Commands");
		if (ImGui::Button("Lloyd", ImVec2(-1, 0))) {
			cellogram::lloyd_relaxation(scene.detected, scene.fixed, 1);//, &scene.border_mesh);
			update_mesh(scene.detected, mesh_);
			mesh_gfx_.set_mesh(&mesh_);
		}
	}

	virtual void draw_scene() override {
		SimpleMeshApplication::draw_scene();
		if (show_border_) {
			draw_border();
		}
	}

	void draw_border() {
		glupSetColor3f(GLUP_FRONT_AND_BACK_COLOR, 0.0, 0.0, 0.0);
		glupSetMeshWidth(4);
		glupBegin(GLUP_LINES);
		for(GEO::index_t i=0; i < scene.border.size(); ++i) {
			glupVertex(scene.border[i]);
			glupVertex(scene.border[(i+1)%scene.border.size()]);
		}
		glupEnd();
	}

	void draw_reference() {

	}

};

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
	setenv("GEO_NO_SIGNAL_HANDLERS", "1", 1);
	GEO::initialize();
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::import_arg_group("algo");
	GEO::CmdLine::import_arg_group("gfx");
	CellogramApplication app(argc, argv, "<filename>");
	GEO::CmdLine::set_arg("gfx:geometry", "1024x1024");
	app.start();
	return 0;
}
