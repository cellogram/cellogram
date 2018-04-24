project(delaunay_test)

add_executable(${PROJECT_NAME}
	delaunay.cpp
)

add_executable(PointsDetection
	pts_detection.cpp
)

add_executable(Untangler
	pts_untangler.cpp
)

add_executable(tcdf_test
	tcdf_test.cpp
)

# Dependencies
target_link_libraries(${PROJECT_NAME}
	tinyfiledialogs
	geogram
	geogram_gfx
	gurobi_solver
)

target_link_libraries(PointsDetection
	cellogram
)

target_link_libraries(Untangler
	cellogram PointsUntangler
)

target_link_libraries(tcdf_test
	tcdf
)

# Enable extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})
target_compile_options(PointsDetection PRIVATE ${ALL_WARNINGS})
target_compile_options(Untangler PRIVATE ${ALL_WARNINGS})
target_compile_options(tcdf_test PRIVATE ${ALL_WARNINGS})

# Output directory for binary
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set_target_properties(PointsDetection PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set_target_properties(Untangler PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set_target_properties(tcdf_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# Data folder
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../data/")
target_compile_definitions(${PROJECT_NAME} PRIVATE -DDATA_DIR=\"${DATA_DIR}\")
target_compile_definitions(PointsDetection PRIVATE -DDATA_DIR=\"${DATA_DIR}\")
target_compile_definitions(Untangler PRIVATE -DDATA_DIR=\"${DATA_DIR}\")
target_compile_definitions(tcdf_test PRIVATE -DDATA_DIR=\"${DATA_DIR}\")
