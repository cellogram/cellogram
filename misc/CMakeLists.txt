cmake_minimum_required(VERSION 3.1)

################################################################################

# Data folder
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../data/")

function(add_application APP_SOURCE)
	# Add executable
	get_filename_component(APP_NAME ${APP_SOURCE} NAME_WE)
	add_executable(${APP_NAME} ${APP_SOURCE})
	message(STATUS "Compiling single-source application: ${APP_NAME}")

	# Dependencies
	target_link_libraries(${APP_NAME} PRIVATE ${ARGN})

	# Output directory for binaries
	set_target_properties(${APP_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

	# Enable extra warnings
	target_compile_options(${APP_NAME} PRIVATE ${ALL_WARNINGS})

	# Output directory for binary
	set_target_properties(${APP_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

	# Data folder
	target_compile_definitions(${APP_NAME} PRIVATE -DDATA_DIR=\"${DATA_DIR}\")

	set_target_properties(${APP_NAME} PROPERTIES FOLDER test)
endfunction()

################################################################################

add_application(pts_detection.cpp cellogram)
add_application(tcdf_test.cpp tcdf)
add_application(tif_test.cpp cellogram)
add_application(image_test.cpp cellogram)
add_application(adaptive_test.cpp cellogram)
add_application(pts_untangler.cpp cellogram)
