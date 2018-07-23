////////////////////////////////////////////////////////////////////////////////
#include "StringUtils.h"
#include <iomanip>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {
namespace StringUtils {

// -----------------------------------------------------------------------------


	int cellogram_mkdir(const std::string &path) {
		int nError;
#if defined(_WIN32)
		std::wstring widestr = std::wstring(path.begin(), path.end());
		nError = _wmkdir(widestr.c_str()); // can be used on Windows
#else
		nError = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); // can be used on non-Windows
#endif
		return nError;
	}

// Split a string into tokens
std::vector<std::string> split(const std::string &str, const std::string &delimiters) {
	// Skip delimiters at beginning.
	std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

	std::vector<std::string> tokens;
	while (std::string::npos != pos || std::string::npos != lastPos) {
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}

	return tokens;
}

// Skip comments in a stream
std::istream &skip(std::istream &in, char x) {
	std::string dummy;
	while ((in >> std::ws).peek() ==
		std::char_traits<char>::to_int_type(x))
	{
		std::getline(in, dummy);
	}
	return in;
}

// Tests whether a string starts with a given prefix
bool startswith(const std::string &str, const std::string &prefix) {
	return (str.compare(0, prefix.size(), prefix) == 0);
}

// Tests whether a string ends with a given suffix
bool endswidth(const std::string &str, const std::string &suffix) {
	if (str.length() >= suffix.length()) {
		return (0 == str.compare(str.length() - suffix.length(), suffix.length(), suffix));
	} else {
		return false;
	}
}

// -----------------------------------------------------------------------------

} // namespace StringUtils
} // namespace cellogram
