cmake_minimum_required(VERSION 3.0)

project(phi4-2pt CXX C)

add_library(phi4-common
	action.cpp
	complex.cpp
	metropolis.cpp
	correlators.cpp
	scalar.cpp
	)

find_package(OpenMP)

# Default to "Release" build type.
message(STATUS "Build Type: '${CMAKE_BUILD_TYPE}'")
if(CMAKE_BUILD_TYPE STREQUAL "")
	  message(STATUS "No CMAKE_BUILD_TYPE, assuming release and enabling all optimizations.")
	    set(CMAKE_BUILD_TYPE Release)
    endif()

# Tune for the local architecture because this code will be used where it is
# compiled, in most cases.
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

target_compile_options(phi4-common PUBLIC ${OpenMP_C_FLAGS} --std=c++11)
target_link_libraries(phi4-common PUBLIC ${OpenMP_C_FLAGS} -lm)

add_executable(toytest
	calculate_toytest.cpp
	)
target_link_libraries(toytest PUBLIC phi4-common)

add_executable(corr
	calculate_corr.cpp
	)
target_link_libraries(corr PUBLIC phi4-common)
