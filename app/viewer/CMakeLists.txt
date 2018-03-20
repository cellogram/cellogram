project(viewer)

# Main executable
add_executable(${PROJECT_NAME}
	main.cpp
	UIState.cpp
	UIState.h
)

# Dependencies
target_link_libraries(${PROJECT_NAME} PUBLIC
	igl::opengl_glfw_imgui
	igl::png
	cellogram
	tinyfiledialogs
)

# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})

# Output directory for binary
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
