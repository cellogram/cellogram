////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include "CLI11.hpp"
#include <cellogram/State.h>
#include <cellogram/StringUtils.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
////////////////////////////////////////////////////////////////////////////////

using namespace cellogram;

#ifdef WIN32
int setenv(const char *name, const char *value, int overwrite) {
	int errcode = 0;
	if (!overwrite) {
		size_t envsize = 0;
		errcode = getenv_s(&envsize, NULL, 0, name);
		if (errcode || envsize) return errcode;
	}
	return _putenv_s(name, value);
}
#endif

int main(int argc, char *argv[]) {
	setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);

	GEO::initialize();

	// Import standard command line arguments, and custom ones
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::import_arg_group("pre");
	GEO::CmdLine::import_arg_group("algo");

	// Default arguments
	struct {
		std::string input = DATA_DIR "perfects.png";
		std::string settings = "";
		bool cmd = false;
	} args;

	// Parse arguments
	CLI::App app{"cellogram"};
	app.add_option("input,-i,--input", args.input, "Input image.");
	app.add_option("-s,--settings", args.settings, "Path to json settings");
	app.add_flag("-c,--cmd", args.cmd, "Run without GUI");
	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		return app.exit(e);
	}

	if (args.cmd)
	{
		auto &state = cellogram::State::state();
		state.load_image(args.input);

		if (state.img.size() == 0)
		{
			std::cout << "Image not loaded" << std::endl;
			exit(0);
		}
		state.load_settings(args.settings);

		state.detect_vertices();
		state.untangle();
		state.detect_bad_regions();
		state.resolve_regions();
		state.final_relax();

		if (!state.image_from_pillars)
		{
			state.mesh_3d_uniform();
			state.remesh_3d_adaptive();
		}
		state.analyze_3d_mesh();

		const int index = args.input.find_last_of(".");
		std::string save_dir = args.input.substr(0, index);
		StringUtils::cellogram_mkdir(save_dir);
		state.save(save_dir);
	}
	else
	{
		UIState::ui_state().initialize();
		UIState::ui_state().load_image(args.input);
		UIState::ui_state().state.load_settings(args.settings);
		UIState::ui_state().launch();
	}
	return 0;
}
