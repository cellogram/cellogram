################################################################################

cmake_minimum_required(VERSION 3.3)
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

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	foreach(config ${CMAKE_CONFIGURATION_TYPES})
		string(TOUPPER ${config} config)
		string(REPLACE /MD /MT CMAKE_C_FLAGS_${config} "${CMAKE_C_FLAGS_${config}}")
		string(REPLACE /MD /MT CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
	endforeach()
endif()

################################################################################

# Options for libigl modules
option(LIBIGL_USE_STATIC_LIBRARY          "Use libigl as static library" ON)
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

option(CELLOGRAM_WTIH_UNTANGLER           "Use Marco untangler"          OFF)

################################################################################
# Required libraries
################################################################################

add_subdirectory(${THIRD_PARTY_DIR}/polyfem/cpp)

# Eigen
# Uncomment to use the system's version of Eigen (e.g. to use Eigen 3.3)
# find_package(Eigen3 REQUIRED)
# target_link_libraries(polyfem PUBLIC Eigen3::Eigen)

# libigl
# find_package(LIBIGL REQUIRED)

# Geogram
# include(geogram)

# CLI11
# add_library(cli11 INTERFACE)
# target_include_directories(cli11 SYSTEM INTERFACE ${THIRD_PARTY_DIR}/cli11)

# Clipper
# add_subdirectory(${THIRD_PARTY_DIR}/clipper)

# Json
# add_library(json INTERFACE)
# target_include_directories(json SYSTEM INTERFACE ${THIRD_PARTY_DIR}/json/single_include)

# GSL
add_subdirectory(${THIRD_PARTY_DIR}/gsl)
target_include_directories(gsl PUBLIC ${GSL_BINARY_DIR})

# Tiny file dialogs library
# add_subdirectory(${THIRD_PARTY_DIR}/tinyfiledialogs tinyfiledialogs)

# Tiny format library
# add_library(tinyformat INTERFACE)
# target_include_directories(tinyformat SYSTEM INTERFACE ${THIRD_PARTY_DIR}/tinyformat)

# Gurobi library
find_package(GUROBI)
if(NOT GUROBI_FOUND)
	message(WARNING "Gurobi was not found on this computer, some features will be disabled.")
endif()

# TIFF
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

################################################################################

add_subdirectory(src)
add_subdirectory(app)
add_subdirectory(misc)
add_subdirectory(tests)

################################################################################

set_target_properties(TinyTiff PROPERTIES FOLDER 3rdparty)
set_target_properties(block_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(bspline_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(cblas_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(cdf_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(cheb_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(combination_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(complex_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(const_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(copy-headers PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(deriv_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(dht_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(diff_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(eigen_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(err_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(fft_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(fit_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(glad PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(gsl PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(gslcblas PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(histogram_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(ieee-utils_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(igl_core PROPERTIES FOLDER 3rdparty/igl)
set_target_properties(igl_opengl PROPERTIES FOLDER 3rdparty/igl)
set_target_properties(igl_opengl_glfw PROPERTIES FOLDER 3rdparty/igl)
set_target_properties(igl_opengl_glfw_imgui PROPERTIES FOLDER 3rdparty/igl)
set_target_properties(igl_png PROPERTIES FOLDER 3rdparty/igl)
set_target_properties(igl_stb_image PROPERTIES FOLDER 3rdparty/igl)
set_target_properties(igl_tetgen PROPERTIES FOLDER 3rdparty/igl)
set_target_properties(igl_triangle PROPERTIES FOLDER 3rdparty/igl)
set_target_properties(imgui PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(integration_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(interpolation_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(linalg_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(matrix_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(matrix_test_static PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(min_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(monte_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(multifit_nlinear_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(multifit_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(multilarge_nlinear_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(multilarge_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(multimin_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(multiroots_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(multiset_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(ntuple_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(ode-initval_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(ode-initval2_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(permutation_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(poly_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(polyclipping PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(qrng_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(randist_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(rng_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(roots_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(rstat_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(siman_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(sort_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(spblas_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(specfunc_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(splinalg_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(spmatrix_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(statistics_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(sum_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(sys_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(tetgen PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(tinyfiledialogs PROPERTIES FOLDER 3rdparty)
set_target_properties(triangle PROPERTIES FOLDER 3rdparty)
set_target_properties(vector_test PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(vector_test_static PROPERTIES FOLDER 3rdparty/not_used)
set_target_properties(wavelet_test PROPERTIES FOLDER 3rdparty/not_used)
