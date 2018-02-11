# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( max_cover )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )

if (CGAL_FOUND)
  if (CMAKE_COMPILER_IS_GNUCXX)
    add_definitions(-std=c++11)
  endif()

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  
  include_directories("${PROJECT_BINARY_DIR}")

  add_executable(max_cover main.cpp cover.cpp cover.hpp globals.hpp)

else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()
