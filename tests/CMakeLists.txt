cmake_minimum_required(VERSION 3.1)
################################################################################

set(test_sources
	main.cpp
	InterpTest.cpp
	CylinderTest.cpp
	OptimTest.cpp
)


add_executable(tests ${test_sources})

foreach(source IN ITEMS ${test_sources})
    source_group("tests" FILES "${source}")
endforeach()

target_link_libraries(tests PUBLIC cellogram zebrafish Catch2::Catch2 warnings::all igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui igl::tetgen)

set(DATA_DIR "$CMAKE_CURRENT_SOURCE_DIR/data/")
target_include_directories(tests PUBLIC ${THIRD_PARTY_DIR}/Catch2/single_include/catch2)
#target_compile_definitions(tests PUBLIC -DZEBRA_DATA_DIR=\"${DATA_DIR}\")


# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(tests)
