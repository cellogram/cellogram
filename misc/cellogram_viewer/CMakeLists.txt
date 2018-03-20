project(old_cellogram_viewer)

add_executable(${PROJECT_NAME}
	main.cpp
	file_dialog.cpp
	file_dialog.h
)

# Dependencies
target_link_libraries(${PROJECT_NAME} cellogram tinyfiledialogs geogram geogram_gfx)

# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})

# Output directory for binary
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# Data folder
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../data/")
target_compile_definitions(${PROJECT_NAME} PRIVATE -DDATA_DIR=\"${DATA_DIR}\")
