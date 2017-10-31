################################################################################

add_library(${PROJECT_NAME}
		cellogram/common.cpp
		cellogram/common.h
		cellogram/convex_hull.cpp
		cellogram/convex_hull.h
		cellogram/voronoi.cpp
		cellogram/voronoi.h
)

# Include directory
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# C++11
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)

# Let's get a little bit paranoid
include(SetWarnings)
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Dependencies
target_link_libraries(${PROJECT_NAME} geogram)

################################################################################

add_executable(${PROJECT_NAME}_viewer
		viewer/main.cpp
		viewer/file_dialog.cpp
		viewer/file_dialog.h
)

# Include directory
target_include_directories(${PROJECT_NAME}_viewer PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# C++11
set_target_properties(${PROJECT_NAME}_viewer PROPERTIES CXX_STANDARD 11)
set_target_properties(${PROJECT_NAME}_viewer PROPERTIES CXX_STANDARD_REQUIRED ON)

# Let's get a little bit paranoid
include(SetWarnings)
target_compile_options(${PROJECT_NAME}_viewer PRIVATE ${ALL_WARNINGS})

# Output directory for binary
set_target_properties(${PROJECT_NAME}_viewer
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Data folder
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../data/")
target_compile_definitions(${PROJECT_NAME}_viewer PRIVATE -DDATA_DIR=\"${DATA_DIR}\")

# Dependencies
target_link_libraries(${PROJECT_NAME}_viewer ${PROJECT_NAME} tinyfiledialogs geogram geogram_gfx)

################################################################################
add_subdirectory(gurobiImplementation)