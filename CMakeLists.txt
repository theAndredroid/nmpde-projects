cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(cmake-common.cmake)

add_custom_command(OUTPUT mesh
                  COMMAND mkdir mesh)

add_custom_command(OUTPUT mesh/mesh-square.msh
                  # COMMAND gmsh ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh-square.geo -setnumber h 0.1 -save -o mesh/mesh-square.msh
                  COMMAND gmsh ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh-square.geo -save -o mesh/mesh-square.msh
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh-square.geo
                  DEPENDS mesh)

add_custom_target(mesh_file ALL
                  DEPENDS mesh/mesh-square.msh)

add_executable(exercise-01 src/exercise-01.cpp src/Current.cpp)
deal_ii_setup_target(exercise-01)
