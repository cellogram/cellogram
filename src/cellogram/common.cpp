////////////////////////////////////////////////////////////////////////////////
#include "common.h"
#include <stdexcept>
#include <sstream>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

[[noreturn]] void assertion_failed(
	const std::string& condition_string,
	const std::string& file, int line)
{
	std::ostringstream os;
	os << "Assertion failed: " << condition_string << ".\n";
	os << "File: " << file << ",\n";
	os << "Line: " << line;

	throw std::runtime_error(os.str());
}

} // namespace cellogram
