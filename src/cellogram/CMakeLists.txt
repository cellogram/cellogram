project(cellogram)

add_library(${PROJECT_NAME}
	common.cpp
	common.h
	convex_hull.cpp
	convex_hull.h
	voronoi.cpp
	voronoi.h
)

target_include_directories(${PROJECT_NAME} PUBLIC ..)

# Dependencies
target_link_libraries(${PROJECT_NAME} geogram)

# Enable extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Output directory for binary
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
