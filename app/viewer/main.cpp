////////////////////////////////////////////////////////////////////////////////
#include <CLI/CLI.hpp>
#include "UIState.h"
#include <cellogram/State.h>
#include <cellogram/StringUtils.h>
#include <cellogram/PNGOutput.h>
#include <zebrafish/Logger.hpp>

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
		std::string load_data = "";
		int start_phase = 0;
		int end_phase = 3;
		bool cmd = false;
	} args;

	// Parse arguments
	CLI::App app{"cellogram"};
	app.add_option("input,-i,--input", args.input, "Input image.");
	app.add_option("-s,--settings", args.settings, "Path to json settings");
	app.add_option("-f,--file", args.load_data, "Path to saved data for scene");
	app.add_option("-b,--begin", args.start_phase, "From which phase to run the script");
	app.add_option("-e,--end", args.end_phase, "Until which phase to run the script");
	app.add_flag("-c,--cmd", args.cmd, "Run without GUI");

	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		return app.exit(e);
	}

	if (args.cmd)
	{
		PNGOutput png_output(2);
		auto &state = cellogram::State::state();
		state.load_image(args.input);

		if (state.img.size() == 0)
		{
			std::cout << "Image not loaded" << std::endl;
			exit(0);
		}

		// load data from path
		bool ok = true;
		if (args.load_data.size() > 0)
		{
			std::cout << "Loading data..." << std::endl;
			ok = state.load(args.load_data);
		}
		if (!ok)
		{
			std::cout << "Data not loaded" << std::endl;
			exit(0);
		}
		else
		{
			std::cout << "... data loaded" << std::endl;
		}

		const int index = args.input.find_last_of(".");
		std::string save_dir = args.input.substr(0, index);
		std::cout<<save_dir<<std::endl;
		StringUtils::cellogram_mkdir(save_dir);

		state.load_settings(args.settings);
#ifdef CELLOGRAM_WITH_PNG
		std::string image_extension = ".png";
#else
		std::string image_extension = ".svg";
#endif
		png_output.init(save_dir + "/all" + image_extension, state.img.rows(), state.img.cols());

		png_output.draw_image();

		if(args.end_phase > 0 && args.start_phase < 2){
			std::cout << "Detecting vertices" << std::endl;
			state.detect_vertices();
			png_output.draw_detection();
			state.phase_enumeration = 1;
		}
		if (args.end_phase > 1 && args.start_phase < 3)
		{
			std::cout << "Untangling" << std::endl;
			state.untangle();
			state.detect_bad_regions();
			state.resolve_regions();


			png_output.draw_untangle();

			state.final_relax();
			state.phase_enumeration = 3;
		}
		if (args.end_phase > 2 && args.start_phase < 4) {
			state.phase_enumeration = 4;
		}
		if (args.end_phase > 3 && args.start_phase < 5)
		{
			if (!state.image_from_pillars)
			{
				std::cout << "Meshing volume" << std::endl;
				state.mesh_3d_volume();
				std::cout << "Remeshing adaptively" << std::endl;
				state.remesh_3d_adaptive();
			}
			std::cout << "Running analysis" << std::endl;
			state.analyze_3d_mesh();
			state.phase_enumeration = 5;
		}

		state.save(save_dir);

		png_output.save();
	}  // command line mode
	else
	{
		// logger
		int log_level = 0;
		zebrafish::Logger::init(UIState::ui_state().oss);
		log_level = std::max(0, std::min(6, log_level));
		spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
		spdlog::flush_every(std::chrono::seconds(3));
		zebrafish::logger().info("Cellogram2_gui logger initialized.");

		UIState::ui_state().initialize();
		UIState::ui_state().state.load_settings(args.settings);
		UIState::ui_state().load_image(args.input);
		UIState::ui_state().launch();
	}
	return 0;
}
