################################################################################

cmake_minimum_required(VERSION 3.3)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/StaticRuntime.cmake NO_POLICY_SCOPE)
project(cellogram)

################################################################################

set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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

# C++11/14 features
include(CXXFeatures)

# Export compile flags (used for autocompletion of the C++ code)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# Generate position independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################

# Options for libigl modules
option(LIBIGL_USE_STATIC_LIBRARY          "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR            "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL                   "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO                 "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK                   "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE                 "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM                    "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB                 "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK                  "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL                 "Use OpenGL"                   ON)
option(LIBIGL_WITH_OPENGL_GLFW            "Use GLFW"                     ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI      "Use ImGui"                    ON)
option(LIBIGL_WITH_PNG                    "Use PNG"                      ON)
option(LIBIGL_WITH_PYTHON                 "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN                 "Use Tetgen"                   ON)
option(LIBIGL_WITH_TRIANGLE               "Use Triangle"                 ON)
option(LIBIGL_WITH_VIEWER                 "Use OpenGL viewer"            ON)
option(LIBIGL_WITH_XML                    "Use XML"                      OFF)
option(CELLOGRAM_WITH_PNG           	  "Use PNG output"				 ON)
option(CELLOGRAM_WITH_MISC           	  "Use build misc targets"		 OFF)

################################################################################
# Download external dependencies
################################################################################

include(CellogramDownloadExternal)

################################################################################
# Required libraries
################################################################################

cellogram_download_polyfem()
add_subdirectory(${THIRD_PARTY_DIR}/polyfem)

if(CELLOGRAM_WITH_PNG)
	find_package(PNG 1.2.9)
	find_package(ZLIB)
	if(NOT ${ZLIB_FOUND} OR NOT ${PNG_FOUND})
		set(CELLOGRAM_WITH_PNG OFF)
		message(WARNING "disabling png output, coudn't find png or zlip")
	endif()
endif()

if(CELLOGRAM_WITH_PNG)
	cellogram_download_pngwriter()
	add_subdirectory(${THIRD_PARTY_DIR}/pngwriter)
endif()

# Eigen
# Uncomment to use the system's version of Eigen (e.g. to use Eigen 3.3)
# find_package(Eigen3 REQUIRED)
# target_link_libraries(polyfem PUBLIC Eigen3::Eigen)

# GSL
cellogram_download_gsl()
SET(GSL_DISABLE_WARNINGS ON)
SET(GSL_DISABLE_TESTS ON)
add_subdirectory(${THIRD_PARTY_DIR}/gsl)
target_include_directories(gsl SYSTEM PUBLIC ${GSL_BINARY_DIR})

# MMG
cellogram_download_mmg()
option(BUILD_TESTING "Enable/Disable continuous integration" OFF)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
add_subdirectory(${THIRD_PARTY_DIR}/mmg)
add_library(mmg::mmg ALIAS libmmg_a)
add_library(mmg::mmgs ALIAS libmmgs_a)
add_library(mmg::mmg2d ALIAS libmmg2d_a)
add_library(mmg::mmg3d ALIAS libmmg3d_a)

# TIFF
cellogram_download_tinytiff()
add_library(TinyTiff
	${THIRD_PARTY_DIR}/TinyTIFF/tinytiffhighrestimer.cpp
	${THIRD_PARTY_DIR}/TinyTIFF/tinytiffhighrestimer.h
	${THIRD_PARTY_DIR}/TinyTIFF/tinytiffreader.cpp
	${THIRD_PARTY_DIR}/TinyTIFF/tinytiffreader.h
	${THIRD_PARTY_DIR}/TinyTIFF/tinytiffwriter.cpp
	${THIRD_PARTY_DIR}/TinyTIFF/tinytiffwriter.h
	# ${THIRD_PARTY_DIR}/TinyTIFF/libtiff_tools/libtiff_tools.cpp
	# ${THIRD_PARTY_DIR}/TinyTIFF/libtiff_tools/libtiff_tools.h
)
target_include_directories(TinyTiff SYSTEM INTERFACE ${THIRD_PARTY_DIR}/TinyTIFF)


if(NOT TARGET tinyfiledialogs::tinyfiledialogs)
	cellogram_download_tinyfiledialogs()
	set(TINYFILEDIALOGS_DIR "${THIRD_PARTY_DIR}/tinyfiledialogs")
	add_library(tinyfiledialogs_tinyfiledialogs ${TINYFILEDIALOGS_DIR}/tinyfiledialogs.c)
	add_library(tinyfiledialogs::tinyfiledialogs ALIAS tinyfiledialogs_tinyfiledialogs)
	target_include_directories(tinyfiledialogs_tinyfiledialogs SYSTEM INTERFACE ${TINYFILEDIALOGS_DIR})
	set_target_properties(tinyfiledialogs_tinyfiledialogs PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

################################################################################

add_subdirectory(src)


if(NOT TARGET CLI11::CLI11)
    cellogram_download_cli11()
    add_subdirectory(${THIRD_PARTY_DIR}/cli11)
endif()
add_subdirectory(app)
if(CELLOGRAM_WITH_MISC)
	add_subdirectory(misc)
endif()

################################################################################

set_target_properties(TinyTiff PROPERTIES FOLDER 3rdparty)
set_target_properties(copy-headers PROPERTIES FOLDER 3rdparty/not_used)
# set_target_properties(glad PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(gsl PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(gslcblas PROPERTIES FOLDER 3rdparty/not_used)
# set_target_properties(igl PROPERTIES FOLDER 3rdparty/igl)
# set_target_properties(igl_opengl PROPERTIES FOLDER 3rdparty/igl)
# set_target_properties(igl_opengl_glfw PROPERTIES FOLDER 3rdparty/igl)
# set_target_properties(igl_opengl_glfw_imgui PROPERTIES FOLDER 3rdparty/igl)
# set_target_properties(igl_png PROPERTIES FOLDER 3rdparty/igl)
# set_target_properties(igl_stb_image PROPERTIES FOLDER 3rdparty/igl)
# set_target_properties(igl_tetgen PROPERTIES FOLDER 3rdparty/igl)
# set_target_properties(igl_triangle PROPERTIES FOLDER 3rdparty/igl)
# set_target_properties(imgui PROPERTIES FOLDER 3rdparty/not_used)
# set_target_properties(polyclipping PROPERTIES FOLDER 3rdparty/not_used)
# set_target_properties(tetgen PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(tinyfiledialogs_tinyfiledialogs PROPERTIES FOLDER 3rdparty)
# set_target_properties(triangle PROPERTIES FOLDER 3rdparty)
