if(CELLOGRAM_WTIH_UNTANGLER)
	project(cellogram_untangler)

	set(source_list
		grid.cpp
		grid.h
		mesh.cpp
		mesh.h
		mesh_delaunay.cpp
		mesh_io.cpp
		mesh_to_grid.cpp
		my_assert.h
		# orient2D.cpp
		points_untangler.cpp
		points_untangler.h
		vec2.h
	)

	add_library(${PROJECT_NAME} ${source_list})
	target_include_directories(${PROJECT_NAME} PUBLIC ..)
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DCELLOGRAM_WITH_UNTANGLER)
	target_link_libraries(${PROJECT_NAME} igl::core)

	# Extra warnings
	target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

	# Use C++14
	target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})

	foreach(source IN ITEMS ${source_list})
	    get_filename_component(source_path "${source}" PATH)
	    string(REPLACE "/" "\\" group_path "${source_path}")
	    source_group("${group_path}" FILES "${source}")
	endforeach()

	# Alias target
	add_library(cellogram::untangler ALIAS ${PROJECT_NAME})
else()
	# Dummy empty target
	add_library(cellogram_untangler INTERFACE)
	add_library(cellogram::untangler ALIAS cellogram_untangler)
endif()
