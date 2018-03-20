# Main library
add_library(gurobi_solver
	generateQ.cpp
	generateQ.h
	gurobiModel.cpp
	gurobiModel.h
	state.cpp
	state.h
)

target_include_directories(gurobi_solver PUBLIC ..)

target_link_libraries(gurobi_solver igl::core geogram gurobi)

# Example test
add_executable(gurobi_example
	main.cpp
)

target_link_libraries(gurobi_example gurobi_solver)

# Data folder for the example
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/")
target_compile_definitions(gurobi_example PUBLIC -DDATA_DIR=\"${DATA_DIR}\")
