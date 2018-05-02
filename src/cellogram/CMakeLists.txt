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
	region_grow.cpp
	region_grow.h
	laplace_energy.cpp
	laplace_energy.h
	load_points.cpp
	load_points.h
	Mesh.cpp
	Mesh.h
	MeshUtils.cpp
	MeshUtils.h
	navigation.cpp
	navigation.h
	PolygonUtils.cpp
	PolygonUtils.h
	State.cpp
	State.h
	StringUtils.cpp
	StringUtils.h
	point_source_detection.cpp
	point_source_detection.h
	Region.cpp
	Region.h
	tri2hex.cpp
	tri2hex.h
	vertex.cpp
	vertex.h
	vertex_degree.cpp
	vertex_degree.h
	voronoi.cpp
	voronoi.h
	fitGaussian2D.cpp
	fitGaussian2D.h
	tif_reader.cpp
	tif_reader.h
)

add_library(${PROJECT_NAME} ${source_list})

target_include_directories(${PROJECT_NAME} PUBLIC ..)

# Dependencies
target_link_libraries(${PROJECT_NAME}
	geogram
	igl::core
	igl::triangle
	igl::opengl_glfw_imgui
	json
	polyclipping
	tinyformat
	gurobi_solver
	tcdf
	points_untangler
	gsl
	TinyTiff
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