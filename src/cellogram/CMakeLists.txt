project(cellogram)

add_library(${PROJECT_NAME}
	common.cpp
	common.h
	convex_hull.cpp
	convex_hull.h
	load_points.cpp
	load_points.h
	MeshUtils.cpp
	MeshUtils.h
	State.cpp
	State.h
	StringUtils.cpp
	StringUtils.h
	voronoi.cpp
	voronoi.h
)

target_include_directories(${PROJECT_NAME} PUBLIC ..)

# Dependencies
target_link_libraries(${PROJECT_NAME}
	geogram
	gurobi_solver
	igl::core
	json
	polyclipping
	tinyformat
)

# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})

