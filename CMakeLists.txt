cmake_minimum_required(VERSION 2.6.2)

cmake_policy(SET CMP0043 NEW)
cmake_policy(SET CMP0072 NEW)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")

project(MoA)

set(CMAKE_CXX_STANDARD 17)

###
## WX
#
if(WIN32)
	set(wxWidgets_ROOT_DIR ${CMAKE_SOURCE_DIR}/VC++/wxWidgets)
	set(wxWidgets_LIB_DIR ${CMAKE_SOURCE_DIR}/VC++/wxWidgets/lib/vc_lib)
	set(wxWidgets_CONFIGURATION msw)
endif(WIN32)

FIND_PACKAGE(wxWidgets REQUIRED base core gl)
INCLUDE(${wxWidgets_USE_FILE})

###
## GLUT
#
if(WIN32)
	set(GLUT_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/VC++/freeglut/include)
	set(GLUT_LIBRARIES ${CMAKE_SOURCE_DIR}/VC++/freeglut/lib/freeglut_static.lib ${CMAKE_SOURCE_DIR}/VC++/freeglut/lib/freeglut_staticd.lib)
	include_directories(${GLUT_INCLUDE_DIR})
else(WIN32)
	FIND_PACKAGE(GLUT REQUIRED)
endif(WIN32)

###
## GLEW
#
if(WIN32)
	set(GLEW_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/VC++/glew/include)
	set(GLEW_LIBRARIES ${CMAKE_SOURCE_DIR}/VC++/glew/lib/Release/x64/glew32.lib)
	include_directories(${GLEW_INCLUDE_DIRS})
else(WIN32)
	FIND_PACKAGE(GLEW REQUIRED)
endif(WIN32)

FIND_PACKAGE(OpenGL REQUIRED)

###
## Output paths for the executables and libraries
#
if(NOT EXECUTABLE_OUTPUT_PATH)
  set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE PATH
    "Single output directory for building all executables.")
endif(NOT EXECUTABLE_OUTPUT_PATH)

if(NOT LIBRARY_OUTPUT_PATH)
  set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib CACHE PATH
    "Single output directory for building all libraries.")
endif(NOT LIBRARY_OUTPUT_PATH)

mark_as_advanced(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH CMAKE_INSTALL_PREFIX)

set(MOA_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/include)
set(MOA_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(MOA_BINARY_DIR ${CMAKE_BINARY_DIR})
set(MOA_LIBRARY_DIR ${LIBRARY_OUTPUT_PATH})
set(MOA_EXECUTABLE_DIR ${EXECUTABLE_OUTPUT_PATH})

###
## Boolean values to build the different projects
#
option(BUILD_LAB1 "Build Lab 1 - Mesh" ON)
option(BUILD_LAB2 "Build Lab 2 - Decimation" ON)
option(BUILD_LAB3 "Build Lab 3 - Subdivision" ON)
option(BUILD_LAB4 "Build Lab 4 - Implicits" ON)
option(BUILD_LAB5 "Build Lab 5 - Levelsets" ON)
option(BUILD_LAB6 "Build Lab 6 - Fluids" ON)

if(BUILD_LAB6)
	set(BUILD_LAB1 ON)
	set(BUILD_LAB2 ON)
	set(BUILD_LAB3 ON)
	set(BUILD_LAB4 ON)
	set(BUILD_LAB5 ON)
endif(BUILD_LAB6)

if(BUILD_LAB5)
	set(BUILD_LAB1 ON)
	set(BUILD_LAB2 ON)
	set(BUILD_LAB3 ON)
	set(BUILD_LAB4 ON)
endif(BUILD_LAB5)

if(BUILD_LAB4)
	set(BUILD_LAB1 ON)
	set(BUILD_LAB2 ON)
	set(BUILD_LAB3 ON)
endif(BUILD_LAB4)

if(BUILD_LAB3)
	set(BUILD_LAB1 ON)
	set(BUILD_LAB2 ON)
endif(BUILD_LAB3)

if(BUILD_LAB2)
	set(BUILD_LAB1 ON)
endif(BUILD_LAB2)

###
## Add pre-processor definitions for each lab
#
if(BUILD_LAB1)
	add_definitions(-DLAB1)
endif(BUILD_LAB1)

if(BUILD_LAB2)
	add_definitions(-DLAB2)
endif(BUILD_LAB2)

if(BUILD_LAB3)
	add_definitions(-DLAB3)
endif(BUILD_LAB3)

if(BUILD_LAB4)
	add_definitions(-DLAB4)
endif(BUILD_LAB4)

if(BUILD_LAB5)
	add_definitions(-DLAB5)
endif(BUILD_LAB5)

if(BUILD_LAB6)
	add_definitions(-DLAB6)
endif(BUILD_LAB6)

###
## Build shared libs or static libs
#
option(BUILD_SHARED_LIBS "Build shared libs or static libs" ON)
mark_as_advanced(BUILD_SHARED_LIBS)

###
## Build modules
#
set(DECIMATION )
set(FLUID )
set(GEOMETRY )
set(GUI )
set(LEVELSET )
set(MY_MATH )
set(SUBDIVISION )
set(UTIL )

include(Decimation/CMakeLists.txt)
include(Fluid/CMakeLists.txt)
include(Geometry/CMakeLists.txt)
include(GUI/CMakeLists.txt)
include(Levelset/CMakeLists.txt)
include(Math/CMakeLists.txt)
include(Subdivision/CMakeLists.txt)
include(Util/CMakeLists.txt)

set(SOURCE
	${DECIMATION}
	${FLUID}
	${GEOMETRY}
	${GUI}
	${LEVELSET}
	${MY_MATH}
	${SUBDIVISION}
	${UTIL}
	main.cpp
)

add_definitions(-DFREEGLUT_STATIC)
add_definitions(-D_USE_MATH_DEFINES)
add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions(-D_SCL_SECURE_NO_WARNINGS)

if(WIN32)
	set(CMAKE_EXE_LINKER_FLAGS /NODEFAULTLIB:LIBCMT)
	set(CMAKE_EXE_LINKER_FLAGS /NODEFAULTLIB:MSVCRT)
endif(WIN32)

source_group("Decimation" FILES ${DECIMATION})
source_group("Fluid" FILES ${FLUID})
source_group("Geometry" FILES ${GEOMETRY})
source_group("GUI" FILES ${GUI} main.cpp)
source_group("Levelset" FILES ${LEVELSET})
source_group("Math" FILES ${MY_MATH})
source_group("Subdivision" FILES ${SUBDIVISION})
source_group("Util" FILES ${UTIL})

include_directories(
	.
	./
	./VC++/glm/
)

###
## APPLICATION
#
if(WIN32)
	add_executable(MoA WIN32 ${SOURCE})
elseif(APPLE)
	add_executable(MoA MACOSX_BUNDLE ${SOURCE})
else(WIN32)
	add_executable(MoA ${SOURCE})
endif(WIN32)

if(WIN32)
    if(NOT ${CMAKE_VERSION} VERSION_LESS "3.6")
        set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT MoA)
    endif()
endif()

TARGET_LINK_LIBRARIES(MoA ${wxWidgets_LIBRARIES})
TARGET_LINK_LIBRARIES(MoA ${GLUT_LIBRARIES})
TARGET_LINK_LIBRARIES(MoA ${OPENGL_LIBRARIES})

if(WIN32)
	TARGET_LINK_LIBRARIES(MoA optimized msvcrt.lib)
	TARGET_LINK_LIBRARIES(MoA optimized msvcmrt.lib)
endif()

if(NOT WIN32)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated -Wno-deprecated-declarations -Wno-unused-result")
endif()
