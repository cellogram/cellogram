project(tcdf)

set(source_list
	ctcdf.cpp
	ctcdf.h
	ctcdf_emxAPI.cpp
	ctcdf_emxAPI.h
	ctcdf_emxutil.cpp
	ctcdf_emxutil.h
	ctcdf_initialize.cpp
	ctcdf_initialize.h
	ctcdf_terminate.cpp
	ctcdf_terminate.h
	ctcdf_types.h
	gammaln.cpp
	gammaln.h
	tcdf.cpp
	tcdf.h
	rtGetInf.cpp
	rtGetInf.h
	rtGetNaN.cpp
	rtGetNaN.h
	rt_nonfinite.cpp
	rt_nonfinite.h
	rtwtypes.h
)

add_library(${PROJECT_NAME} ${source_list})

# Dependencies
target_link_libraries(${PROJECT_NAME}
	igl::core
)

target_include_directories(${PROJECT_NAME} PUBLIC ..)


# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})


foreach(source IN ITEMS ${source_list})
    get_filename_component(source_path "${source}" PATH)
    string(REPLACE "/" "\\" group_path "${source_path}")
    source_group("${group_path}" FILES "${source}")
endforeach()