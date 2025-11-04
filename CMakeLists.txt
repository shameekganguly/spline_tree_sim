
cmake_minimum_required (VERSION 3.10.0)
project (TREE_SIM)
SET(CMAKE_CXX_FLAGS "-std=c++14 -I/usr/include -I/usr/local/include -Wno-deprecated-declarations")

# set common source
set (COMMON_SOURCE
	${PROJECT_SOURCE_DIR}/src/QuadraticSplineVisual.cpp
	${PROJECT_SOURCE_DIR}/src/QuadraticSplineKinematic.cpp
	${PROJECT_SOURCE_DIR}/src/TreeKinematic.cpp
	${PROJECT_SOURCE_DIR}/src/TreeDynamics.cc
	${PROJECT_SOURCE_DIR}/src/TreeParser.cc
	${PROJECT_SOURCE_DIR}/src/TreeVisual.cc
	${PROJECT_SOURCE_DIR}/src/LCPSolver.cpp
	${PROJECT_SOURCE_DIR}/src/collision/Cursor.cc
	${PROJECT_SOURCE_DIR}/src/collision/Sphere.cc
	${PROJECT_SOURCE_DIR}/src/collision/SplineContact.cc
	${PROJECT_SOURCE_DIR}/src/graphics/Graphics.cc
	# ${PROJECT_SOURCE_DIR}/src/timer/LoopTimer.h
	# ${PROJECT_SOURCE_DIR}/src/force_sensor/ForceSensorSim.h
	# ${PROJECT_SOURCE_DIR}/src/force_sensor/ForceSensorDisplay.h
	)
include_directories(${PROJECT_SOURCE_DIR}/src)

# - CHAI3D
find_package(CHAI3D REQUIRED)
include_directories(${CHAI3D_INCLUDE_DIRS})
add_definitions(${CHAI3D_DEFINITIONS})
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
  add_definitions(-D_MACOSX -DMACOSX)
elseif (CMAKE_SYSTEM_NAME MATCHES Linux)
  add_definitions(-D_LINUX -DLINUX)
endif () #NOTE: Windows is currently unsupported
link_directories(${CHAI3D_LIBRARY_DIRS})

# - glfw3
find_package(glfw3 QUIET)
find_library(GLFW_LIBRARY glfw)

# - tinyxml
find_library(TINYXML_LIBRARIES tinyxml2)

set(COMMON_LIBRARIES
	${CHAI3D_LIBRARIES}
	${GLFW_LIBRARY}
	${TINYXML_LIBRARIES}
	)

# add apps
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY                ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG          ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL     ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE        ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO ${PROJECT_SOURCE_DIR}/bin)
set (BINARY_DIR ${PROJECT_SOURCE_DIR}/bin)

# add_subdirectory(01-branch)
# add_subdirectory(02-simple_tree)
# add_subdirectory(03-bare_tree)
add_subdirectory(04-cable)
