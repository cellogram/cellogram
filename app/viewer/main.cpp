////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include "CLI11.hpp"
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
////////////////////////////////////////////////////////////////////////////////

using namespace cellogram;

int main(int argc, char *argv[]) {
#ifndef WIN32
	setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

	GEO::initialize();

	// Import standard command line arguments, and custom ones
	GEO::CmdLine::import_arg_group("standard");
	GEO::CmdLine::import_arg_group("pre");
	GEO::CmdLine::import_arg_group("algo");

	// Default arguments
	struct {
		std::string input = DATA_DIR "small2.xyz";
	} args;

	// Parse arguments
	CLI::App app{"cellogram"};
	app.add_option("input,-i,--input", args.input, "Output points.");
	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		return app.exit(e);
	}

	UIState::ui_state().initialize();
	UIState::ui_state().load(args.input);
	UIState::ui_state().launch();

	return 0;
}
