project(cellogram)

add_library(${PROJECT_NAME}
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
	mesh_solver.cpp
	mesh_solver.h
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
	tri2hex.cpp
	tri2hex.h
	vertex.cpp
	vertex.h
	vertex_degree.cpp
	vertex_degree.h
	voronoi.cpp
	voronoi.h
)

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
)

# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})

