#pragma once

////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <vector>
#include <string>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

namespace StringUtils {

	// Split a string into tokens
	std::vector<std::string> split(const std::string &str, const std::string &delimiters = " ");

	// Skip comments in a stream
	std::istream &skip(std::istream &in, char x = '#');

	// Tests whether a string starts with a given prefix
	bool startswith(const std::string &str, const std::string &prefix);

	// Tests whether a string ends with a given suffix
	bool endswidth(const std::string &str, const std::string &suffix);

}

} // namespace cellogram
