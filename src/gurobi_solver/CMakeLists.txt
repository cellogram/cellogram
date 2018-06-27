if(TARGET gurobi::gurobi)
	# Main library
	add_library(cellogram_gurobi
		generateQ.cpp
		generateQ.h
		gurobiModel.cpp
		gurobiModel.h
		state.cpp
		state.h
	)

	target_include_directories(cellogram_gurobi PUBLIC ..)
	target_compile_definitions(cellogram_gurobi PUBLIC -DCELLOGRAM_WITH_GUROBI)
	target_link_libraries(cellogram_gurobi igl::core geogram gurobi::gurobi)

	# Example test
	add_executable(gurobi_example
		main.cpp
	)

	target_link_libraries(gurobi_example cellogram_gurobi)

	# Data folder for the example
	set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/")
	target_compile_definitions(gurobi_example PUBLIC -DDATA_DIR=\"${DATA_DIR}\")

	# Alias target
	add_library(cellogram::gurobi ALIAS cellogram_gurobi)
else()
	# Dummy empty target
	add_library(cellogram_gurobi INTERFACE)
	add_library(cellogram::gurobi ALIAS cellogram_gurobi)
endif()
