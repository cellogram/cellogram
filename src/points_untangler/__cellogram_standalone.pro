TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

EIGENPATH = "C:/projects/libraries/libigl"

INCLUDEPATH += $$EIGENPATH/include
INCLUDEPATH += $$EIGENPATH/external/eigen

# DEFINES += IGL_STATIC_LIBRARY

# SOURCES += main.cpp
SOURCES +=\
    mesh.cpp \
    grid.cpp \
    mesh_to_grid.cpp \
    mesh_delaunay.cpp \
    __main.cpp \
    points_untangler.cpp \
    mesh_io.cpp


#  $$EIGENPATH/include/igl/orient2D.cpp \


HEADERS += \
    mesh.h \
    vec2.h \
    grid.h \
    my_assert.h


