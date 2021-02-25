project(zebrafish)
set(source_list
	autodiff.h
	Bspline.cpp
	Bspline.h
	Common.h
	Cylinder.cpp
	Cylinder.h
	Logger.cpp
	Logger.hpp
	Quad.ipp
	Quadrature.cpp
	Quadrature.h
	Quantile.cpp
	Quantile.h
)
add_library(${PROJECT_NAME} ${source_list})

target_include_directories(${PROJECT_NAME} PUBLIC ..)

# Dependencies
target_link_libraries(${PROJECT_NAME} PUBLIC
	igl::core
	polyfem
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

