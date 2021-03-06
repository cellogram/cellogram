project(points_untangler)

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

# Dependencies
target_link_libraries(${PROJECT_NAME}
	igl::core
)

# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})


foreach(source IN ITEMS ${source_list})
    get_filename_component(source_path "${source}" PATH)
    string(REPLACE "/" "\\" group_path "${source_path}")
    source_group("${group_path}" FILES "${source}")
endforeach()


add_library(cellogram::untangler ALIAS ${PROJECT_NAME})