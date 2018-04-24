if(TARGET gurobi)
	add_subdirectory(gurobi_solver)
endif()

add_subdirectory(pointsUntangler)
add_subdirectory(tcdf)

add_subdirectory(cellogram)