project(delaunay_test)

add_executable(${PROJECT_NAME}
	delaunay.cpp
)

# Dependencies
target_link_libraries(${PROJECT_NAME}
	tinyfiledialogs
	geogram
	geogram_gfx
	gurobi_solver
)

# Enable extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Output directory for binary
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# Data folder
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../data/")
target_compile_definitions(${PROJECT_NAME} PRIVATE -DDATA_DIR=\"${DATA_DIR}\")
