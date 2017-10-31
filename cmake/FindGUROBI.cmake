# - Try to find GUROBI
# Once done this will define
#
#  GUROBI_FOUND - System has Gurobi
#  gurobi - Interface target to link against Gurobi

if(GUROBI_FOUND)
	return()
endif()

find_path(GUROBI_INCLUDE_DIR
		NAMES gurobi_c++.h
		PATHS "$ENV{GUROBI_HOME}/include"
				"/Library/gurobi502/mac64/include"
				"C:\\libs\\gurobi502\\include"
				"C:\\gurobi751\\win64\\include"
)

find_library(GUROBI_LIBRARY
			NAMES gurobi
					gurobi45
					gurobi46
					gurobi50
					gurobi51
					gurobi55
					gurobi56
					gurobi60
					gurobi75
			PATHS "$ENV{GUROBI_HOME}/lib"
					"/Library/gurobi502/mac64/lib"
					"C:\\libs\\gurobi502\\lib"
					"C:\\gurobi751\\win64\\lib"
)

find_library(GUROBI_CXX_LIBRARY
			NAMES gurobi_c++
					gurobi_c++mtd2017
					gurobi_c++mtd2015
			PATHS "$ENV{GUROBI_HOME}/lib"
					"/Library/gurobi502/mac64/lib"
					"C:\\libs\\gurobi502\\lib"
					"C:\\gurobi751\\win64\\lib"
)

set(GUROBI_INCLUDE_DIRS "${GUROBI_INCLUDE_DIR}")
set(GUROBI_LIBRARIES "${GUROBI_LIBRARY};${GUROBI_CXX_LIBRARY}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG
	GUROBI_LIBRARY GUROBI_CXX_LIBRARY GUROBI_INCLUDE_DIR)

mark_as_advanced(GUROBI_INCLUDE_DIR GUROBI_LIBRARY GUROBI_CXX_LIBRARY)

# Interface target
add_library(gurobi INTERFACE)
target_include_directories(gurobi INTERFACE ${GUROBI_INCLUDE_DIRS})
target_link_libraries(gurobi INTERFACE ${GUROBI_LIBRARIES})