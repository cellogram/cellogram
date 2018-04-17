////////////////////////////////////////////////////////////////////////////////
#include "UIState.h"
#include "CLI11.hpp"
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
		std::string input = DATA_DIR "2_points.xyz";
		std::string param = DATA_DIR "2_param.txt";
	} args;

	// Parse arguments
	CLI::App app{"cellogram"};
	app.add_option("input,-i,--input", args.input, "Output points.");
	app.add_option("param,-p,--param", args.param, "Output params.");
	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		return app.exit(e);
	}

	UIState::ui_state().initialize();
	UIState::ui_state().load(args.input);
	UIState::ui_state().load_param(args.param);
	UIState::ui_state().launch();

	return 0;
}
