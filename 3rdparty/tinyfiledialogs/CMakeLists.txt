cmake_minimum_required(VERSION 3.0)
project(tinyfiledialogs)

if(NOT TARGET tinyfiledialogs)
	add_library(tinyfiledialogs tinyfiledialogs.c)
	target_include_directories(tinyfiledialogs SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
	set_target_properties(tinyfiledialogs PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
