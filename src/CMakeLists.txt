if(TARGET gurobi)
	add_subdirectory(gurobi_solver)
endif()

if(CELLOGRAM_WTIH_UNTANGLER)
	add_subdirectory(pointsUntangler)
endif()
add_subdirectory(tcdf)

add_subdirectory(cellogram)