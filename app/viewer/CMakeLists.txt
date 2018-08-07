project(viewer)

# Resource file (icons etc.)
set(EXTRA_SOURCE "")
if (APPLE)
	set(BUILD_TYPE MACOSX_BUNDLE)
	set(EXTRA_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/resources/im.icns")
elseif(WIN32)
	set(EXTRA_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/resources/im.rc")
endif()

# Main executable
add_executable(${PROJECT_NAME} ${BUILD_TYPE}
	main.cpp
	FileDialog.cpp
	FileDialog.h
	UIMenu.cpp
	UIState.cpp
	UIState.h

	${EXTRA_SOURCE}
)

# Dependencies
target_link_libraries(${PROJECT_NAME} PUBLIC
	igl::opengl_glfw_imgui
	cellogram
	tinyfiledialogs
	cli11
)

# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})

# Output directory for binary
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# Data folder
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../data/")
target_compile_definitions(${PROJECT_NAME} PRIVATE -DDATA_DIR=\"${DATA_DIR}\")


if (APPLE)
  # Build an application bundle on OSX
  set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME Cellogram)
  set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_GUI_IDENTIFIER com.im.cellogram)
  set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_ICON_FILE im.icns)
  set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/resources/MacOSXBundleInfo.plist.in)
  set_source_files_properties(resources/im.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
endif()







INSTALL(TARGETS ${PROJECT_NAME} BUNDLE DESTINATION . RUNTIME DESTINATION bin)

if(APPLE)
	set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/bundle")
	# set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")

	set(CPACK_GENERATOR "DragNDrop")
	set(CPACK_DMG_FORMAT "UDBZ")
	set(CPACK_DMG_VOLUME_NAME "${PROJECT_NAME}")
	set(CPACK_SYSTEM_NAME "OSX")
	set(CPACK_OSX_PACKAGE_VERSION "10.13")
	set(CPACK_PACKAGE_VENDOR "Us")
	set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
	set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}")
	set(CPACK_PACKAGE_ICON "${EXTRA_SOURCE}")
	set(CPACK_DMG_BACKGROUND_IMAGE "${CMAKE_CURRENT_SOURCE_DIR}/resources/splash.png")
	set(CPACK_DMG_DS_STORE "${CMAKE_CURRENT_SOURCE_DIR}/resources/DS_Store")
	include(CPack)
endif()