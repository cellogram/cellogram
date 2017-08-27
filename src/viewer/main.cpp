////////////////////////////////////////////////////////////////////////////////
#include "viewer/file_dialog.h"
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram_gfx/glup_viewer/glup_viewer_gui.h>
////////////////////////////////////////////////////////////////////////////////

class CellogramApplication : public GEO::SimpleMeshApplication {
public:
	CellogramApplication(int argc, char** argv, const std::string& usage) :
		SimpleMeshApplication(argc, argv, "")
	{
		name_ = "[Float] Cellogram Viewer";
	}

	virtual void draw_load_menu() override {
		if(ImGui::MenuItem("Load...")) {
			std::string filename = FileDialog::openFileName();
		}
	}

	virtual void draw_save_menu() override {
		if(ImGui::MenuItem("Save as...")) {
			std::string filename = FileDialog::saveFileName();
		}
	}

	/**
	 * \brief Draws the application menus.
	 * \details This function overloads
	 *  Application::draw_application_menus(). It can be used to create
	 *  additional menus in the main menu bar.
	 */
	virtual void draw_application_menus() {

	}

};

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv) {
	setenv("GEO_NO_SIGNAL_HANDLERS", "1", 1);
	GEO::initialize();
	CellogramApplication app(argc, argv, "<filename>");
	GEO::CmdLine::set_arg("gfx:geometry", "1024x1024");
	app.start();
	return 0;
}
