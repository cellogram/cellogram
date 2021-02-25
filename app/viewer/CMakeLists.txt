project(CellogramViewer)

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

	IconsFontAwesome5.h
	fa-solid-900.hpp

	${EXTRA_SOURCE}
)

# Dependencies
target_link_libraries(${PROJECT_NAME} PRIVATE
	igl::opengl_glfw_imgui
	cellogram
	zebrafish
	tinyfiledialogs::tinyfiledialogs
	CLI11::CLI11
)

# Extra warnings
target_compile_options(${PROJECT_NAME} PRIVATE ${ALL_WARNINGS})

# Use C++14
target_compile_features(${PROJECT_NAME} PUBLIC ${CXX14_FEATURES})


# Data folder
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../data/")
target_compile_definitions(${PROJECT_NAME} PRIVATE -DDATA_DIR=\"${DATA_DIR}\")


if (APPLE)
	# Build an application bundle on OSX
	set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_BUNDLE_NAME ${PROJECT_NAME})
	set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_EXECUTABLE_NAME ${PROJECT_NAME})
	set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_GUI_IDENTIFIER com.im.cellogram)
	set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_ICON_FILE im.icns)

	set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_SHORT_VERSION_STRING "0.1")
	set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_BUNDLE_VERSION "0.1.0")
	set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_COPYRIGHT "© 2018, Cellogram")


  	set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/resources/MacOSXBundleInfo.plist.in)
	set_source_files_properties(resources/im.icns PROPERTIES MACOSX_PACKAGE_LOCATION Resources)
endif()




# Output directory for binary
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
install(TARGETS ${PROJECT_NAME} DESTINATION .)


# set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/bundle")

if(APPLE)
	set(CPACK_GENERATOR "DragNDrop")
	set(CPACK_DMG_FORMAT "UDBZ")
	set(CPACK_DMG_VOLUME_NAME "${PROJECT_NAME}")
	set(CPACK_BUNDLE_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/resources/MacOSXBundleInfo.plist.in)
	# set(CPACK_OSX_PACKAGE_VERSION "10.13")

	set(CPACK_SYSTEM_NAME "OSX")

	set(CPACK_DMG_DS_STORE "${CMAKE_CURRENT_SOURCE_DIR}/resources/DS_Store")
elseif(WIN32)
	set(CPACK_GENERATOR "NSIS")
    set(CPACK_SOURCE_GENERATOR "ZIP")
else()
	set(CPACK_GENERATOR "DEB")
    set(CPACK_SOURCE_GENERATOR "TGZ")
endif()

set(CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR};${PROJECT_NAME};Unspecified;/")

set(CPACK_PACKAGE_VENDOR "Cellogram")
set(CPACK_PACKAGE_NAME "${PROJECT_NAME}")
set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}")
set(CPACK_PACKAGE_ICON "${EXTRA_SOURCE}")
# set(CPACK_DMG_BACKGROUND_IMAGE "${CMAKE_CURRENT_SOURCE_DIR}/resources/splash.png")

set(CPACK_PACKAGE_ICON "${EXTRA_SOURCE}")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/resources/README.md")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/resources/README.md")
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_CURRENT_SOURCE_DIR}/resources/README.md")

set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/bin/*;${CMAKE_SOURCE_DIR}/lib/*;${CMAKE_SOURCE_DIR}/include/*")
set(CPACK_INSTALLED_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/resources/scripts;scripts")

include(CPack)
