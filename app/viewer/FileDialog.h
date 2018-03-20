#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace FileDialog {

// -----------------------------------------------------------------------------

std::string openFileName(const std::string &defaultPath = "./.*",
	const std::vector<std::string> &filters = {}, const std::string &desc = "");

std::string saveFileName(const std::string &defaultPath = "./.*",
	const std::vector<std::string> &filters = {}, const std::string &desc = "");

// -----------------------------------------------------------------------------

} // namespace FileDialog
