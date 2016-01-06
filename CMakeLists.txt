cmake_minimum_required(VERSION 3.3)
project(Homework5)
set(SOURCE_FILES main.cpp AntTweakBar.h elements.cpp elements.h)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#########################################################
# FIND GLUT
#########################################################
find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
link_directories(${GLUT_LIBRARY_DIRS})
add_definitions(${GLUT_DEFINITIONS})
if(NOT GLUT_FOUND)
    message(ERROR "GLUT not found!")
endif(NOT GLUT_FOUND)
#########################################################
# FIND OPENGL
#########################################################
find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
link_directories(${OpenGL_LIBRARY_DIRS})
add_definitions(${OpenGL_DEFINITIONS})
if(NOT OPENGL_FOUND)
    message(ERROR "OPENGL not found!")
endif(NOT OPENGL_FOUND)
#########################################################
# FIND ANTTWEAKBAR
#########################################################
FIND_PATH(ANT_TWEAK_BAR_INCLUDE_PATH AntTweakBar.h
        PATHS
        ${PROJECT_SOURCE_DIR})

FIND_LIBRARY( ANT_TWEAK_BAR_LIBRARY AntTweakBar
        PATHS
        ${PROJECT_SOURCE_DIR}/lib
        )

SET(ANT_TWEAK_BAR_FOUND "NO")
IF (ANT_TWEAK_BAR_INCLUDE_PATH AND ANT_TWEAK_BAR_LIBRARY)
    SET(ANT_TWEAK_BAR_FOUND "YES")
ENDIF (ANT_TWEAK_BAR_INCLUDE_PATH AND ANT_TWEAK_BAR_LIBRARY)

#INCLUDE(${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(AntTweakBar DEFAULT_MSG ANT_TWEAK_BAR_INCLUDE_PATH ANT_TWEAK_BAR_LIBRARY)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/lib)
#########################################################
# Include Files
#########################################################

add_executable(Homework5 ${SOURCE_FILES})

########################################################
# Linking & stuff
#########################################################

# create the program
target_link_libraries(Homework5 ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${ANT_TWEAK_BAR_LIBRARY})