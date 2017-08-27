################################################################################
# General Information
################################################################################

cmake_minimum_required(VERSION 3.3)
project(cellogram)

################################################################################

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

################################################################################

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# 3rd party libraries
################################################################################

# Geogram library
include(geogram)

# Tiny file dialogs library
add_subdirectory(${THIRD_PARTY_DIR}/tinyfiledialogs tinyfiledialogs)

################################################################################

add_subdirectory(src)
