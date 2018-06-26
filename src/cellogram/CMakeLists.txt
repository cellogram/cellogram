project(cellogram)

set(source_list
	boundary_loop.cpp
	boundary_loop.h
	common.cpp
	common.h
	convex_hull.cpp
	convex_hull.h
	delaunay.cpp
	delaunay.h
	dijkstra.cpp
	dijkstra.h
	fitGaussian2D.cpp
	fitGaussian2D.h
	image_reader.cpp
	image_reader.h
	laplace_energy.cpp
	laplace_energy.h
	load_points.cpp
	load_points.h
	Mesh.cpp
	Mesh.h
	Mesh3d.cpp
	Mesh3d.h
	MeshUtils.cpp
	MeshUtils.h
	navigation.cpp
	navigation.h
	point_source_detection.cpp
	point_source_detection.h
	PolygonUtils.cpp
	PolygonUtils.h
	Region.cpp
	Region.h
	region_grow.cpp
	region_grow.h
	remesh_adaptive.cpp
	remesh_adaptive.h
	State.cpp
	State.h
	StringUtils.cpp
	StringUtils.h
	tri2hex.cpp
	tri2hex.h
	vertex.cpp
	vertex.h
	vertex_degree.cpp
	vertex_degree.h
	voronoi.cpp
	voronoi.h
)

add_library(${PROJECT_NAME} ${source_list})

target_include_directories(${PROJECT_NAME} PUBLIC ..)

# Dependencies
target_link_libraries(${PROJECT_NAME}
	polyfem

	mmg::mmg
	igl::opengl_glfw_imgui #needs to be removed when all ok
	igl::png
	gurobi_solver
	tcdf
	gsl
	TinyTiff
)

if(CELLOGRAM_WTIH_UNTANGLER)
	target_link_libraries(${PROJECT_NAME}
		points_untangler
	)

	target_compile_definitions(${PROJECT_NAME} PUBLIC -DWITH_UNTANGLER)
endif()

# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})


foreach(source IN ITEMS ${source_list})
    get_filename_component(source_path "${source}" PATH)
    string(REPLACE "/" "\\" group_path "${source_path}")
    source_group("${group_path}" FILES "${source}")
endforeach()
