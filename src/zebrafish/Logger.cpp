#include <zebrafish/Logger.hpp>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/ostream_sink.h>
#include <spdlog/details/registry.h>
#include <spdlog/details/thread_pool.h>
#include <memory>
#include <mutex>
#include <iostream>

namespace zebrafish
{
	std::shared_ptr<spdlog::async_logger> Logger::logger_;

	// See https://github.com/gabime/spdlog#asynchronous-logger-with-multi-sinks
	void Logger::init(std::vector<spdlog::sink_ptr> &sinks) {
		auto l = spdlog::get("zebrafish");
		bool had_zebra_logger = l != nullptr;
		if (had_zebra_logger)
			spdlog::drop("zebrafish");

		spdlog::init_thread_pool(8192, 1);
		Logger::logger_ =
			std::make_shared<spdlog::async_logger>(
				"zebrafish",
				sinks.begin(), sinks.end(),
				spdlog::thread_pool(), spdlog::async_overflow_policy::block);
		spdlog::register_logger(logger_);

		if (had_zebra_logger)
			logger().warn("Removed another zebrafish logger");
	}

	void Logger::init(bool use_cout, const std::string &filename, bool truncate) {
		std::vector<spdlog::sink_ptr> sinks;
		if (use_cout) {
			sinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
		}
		if (!filename.empty()) {
			sinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename, truncate));
		}

		init(sinks);
	}


	void Logger::init(std::ostream &os) {
		std::vector<spdlog::sink_ptr> sinks;
		sinks.emplace_back(std::make_shared<spdlog::sinks::ostream_sink_mt>(os, false));

		// file log
		char buffer[80];
		time_t rawTime;
		struct tm *timeInfo;
		time(&rawTime);
		timeInfo = localtime(&rawTime);
		strftime(buffer, 80, "%m_%dT%H_%M", timeInfo);
		std::string filename = "./Zebrafish-log-";
		sinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(filename + buffer + ".log", true));

		// init(sinks);
		auto &registry_inst = spdlog::details::registry::instance();

		// create global thread pool if not already exists..
		std::lock_guard<std::recursive_mutex> tp_lock(registry_inst.tp_mutex());
		auto tp = registry_inst.get_tp();
		if (tp == nullptr) {
			tp = std::make_shared<spdlog::details::thread_pool>(spdlog::details::default_async_q_size, 1);
			registry_inst.set_tp(tp);
		}

		logger_ = std::make_shared<spdlog::async_logger>("zebrafish", sinks.begin(), sinks.end(), std::move(tp), spdlog::async_overflow_policy::block);
		spdlog::register_logger(logger_);
	}

} // namespace zebrafish
