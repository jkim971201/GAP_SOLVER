cmake_minimum_required(VERSION 3.5)

project(GAP_SOLVER LANGUAGES CXX)

find_package(OpenMP REQUIRED)

set(CMAKE_C_COMPILER "g++")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_PREFIX_PATH "/home/tools/ortools")
find_package(ortools REQUIRED CONFIG)

set(PROJECT_SRC 
	src/main.cpp
	src/ILPSolver.cpp
	src/GAPBuilder.cpp
	src/GDSolver.cpp
	src/ADMMSolver.cpp
)

add_executable(${PROJECT_NAME} ${PROJECT_SRC})

target_link_libraries(${PROJECT_NAME} 
	PRIVATE 
	  ortools::ortools
		OpenMP::OpenMP_CXX
)
