################################################################################

add_executable(delaunay
		delaunay.cpp
)

# C++11
set_target_properties(delaunay PROPERTIES CXX_STANDARD 11)
set_target_properties(delaunay PROPERTIES CXX_STANDARD_REQUIRED ON)

# Let's get a little bit paranoid
include(SetWarnings)
target_compile_options(delaunay PRIVATE ${ALL_WARNINGS})

# Output directory for binary
set_target_properties(delaunay
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Data folder
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../data/")
target_compile_definitions(delaunay PRIVATE -DDATA_DIR=\"${DATA_DIR}\")

# Dependencies
target_link_libraries(delaunay tinyfiledialogs geogram geogram_gfx)

################################################################################
