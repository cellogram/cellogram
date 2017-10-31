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

# Build static binaries
set(BUILD_SHARED_LIBS OFF)

# Use folder in Visual Studio
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# Enable extra warnings
include(SetWarnings)

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

# libigl library
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR      "Use AntTweakBar"    OFF)
option(LIBIGL_WITH_CGAL             "Use CGAL"           OFF)
option(LIBIGL_WITH_COMISO           "Use CoMiso"         OFF)
option(LIBIGL_WITH_CORK             "Use Cork"           OFF)
option(LIBIGL_WITH_EMBREE           "Use Embree"         OFF)
option(LIBIGL_WITH_LIM              "Use LIM"            OFF)
option(LIBIGL_WITH_MATLAB           "Use Matlab"         OFF)
option(LIBIGL_WITH_MOSEK            "Use MOSEK"          OFF)
option(LIBIGL_WITH_NANOGUI          "Use Nanogui menu"   OFF)
option(LIBIGL_WITH_OPENGL           "Use OpenGL"         OFF)
option(LIBIGL_WITH_OPENGL_GLFW      "Use GLFW"           OFF)
option(LIBIGL_WITH_PNG              "Use PNG"            OFF)
option(LIBIGL_WITH_PYTHON           "Use Python"         OFF)
option(LIBIGL_WITH_TETGEN           "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE         "Use Triangle"       OFF)
option(LIBIGL_WITH_VIEWER           "Use OpenGL viewer"  OFF)
option(LIBIGL_WITH_XML              "Use XML"            OFF)
find_package(LIBIGL REQUIRED QUIET)

# Gurobi library
find_package(GUROBI REQUIRED)

################################################################################

add_subdirectory(src/cellogram)
add_subdirectory(src/gurobiImplementation)
add_subdirectory(src/delaunayTest)
add_subdirectory(src/cellogramViewer)
