#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <type_traits>
////////////////////////////////////////////////////////////////////////////////

namespace cellogram {

[[noreturn]] void assertion_failed(
	const std::string& condition_string,
	const std::string& file, int line);

} // namespace cellogram

// Shortcut alias
namespace cel = cellogram;

////////////////////////////////////////////////////////////////////////////////

#define cel_assert_msg(x, s) do {                            \
    if (!(x)) {                                              \
        cellogram::assertion_failed((s), __FILE__, __LINE__); \
    }                                                        \
} while(0)

// -----------------------------------------------------------------------------

#define cel_assert(x) do {                                  \
    if(!(x)) {                                              \
        cellogram::assertion_failed(#x, __FILE__, __LINE__); \
    }                                                       \
} while (0)


#ifdef CELLOGRAM_NO_DEBUG

#define cel_debug(x)

#else

#define cel_debug(x) do {                                   \
    if(!(x)) {                                              \
        cellogram::assertion_failed(#x, __FILE__, __LINE__); \
    }                                                       \
} while (0)

#endif

////////////////////////////////////////////////////////////////////////////////
// External libraries
////////////////////////////////////////////////////////////////////////////////

// Json
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// Tinyformat
// #include <tinyformat.h>
