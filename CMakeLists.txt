cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(cmake-common.cmake)

add_custom_command(OUTPUT mesh
                  COMMAND mkdir mesh)

add_custom_command(OUTPUT models
                  COMMAND mkdir models)

set(DEFAULT_BASE_MESH_FILE mesh/mesh-square)

add_custom_command(OUTPUT ${DEFAULT_BASE_MESH_FILE}.geo
                  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/${DEFAULT_BASE_MESH_FILE}.geo ${DEFAULT_BASE_MESH_FILE}.geo
                  COMMAND touch ${DEFAULT_BASE_MESH_FILE}.geo
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${DEFAULT_BASE_MESH_FILE}.geo
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/mesh)

set(MESH_COMPILING_COMMAND gmsh ${DEFAULT_BASE_MESH_FILE}.geo -save)
set(DEFAULT_MESH_SIZE 0.1)

add_custom_command(OUTPUT ${DEFAULT_BASE_MESH_FILE}.msh
                  # COMMAND gmsh ${CMAKE_CURRENT_SOURCE_DIR}/mesh/mesh-square.geo -setnumber h 0.1 -save -o mesh/mesh-square.msh
                  COMMAND ${MESH_COMPILING_COMMAND} -setnumber h ${DEFAULT_MESH_SIZE} -o ${DEFAULT_BASE_MESH_FILE}.msh
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${DEFAULT_BASE_MESH_FILE}.geo)

add_custom_command(OUTPUT models/models.csv
                  COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/models/models.csv models/models.csv 
                  COMMAND touch models/models.csv 
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/models/models.csv
                  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/models)

add_custom_target(mesh_file ALL
                  DEPENDS ${DEFAULT_BASE_MESH_FILE}.msh)


add_custom_target(models_file ALL
                  DEPENDS models/models.csv)

string(JOIN " " MESH_COMPILING_COMMAND_STR ${MESH_COMPILING_COMMAND} "")
add_compile_definitions(MESH_COMPILING_COMMAND="${MESH_COMPILING_COMMAND_STR}"s)
add_compile_definitions(DEFAULT_BASE_MESH_FILE="${DEFAULT_BASE_MESH_FILE}"s)
add_compile_definitions(DEFAULT_MESH_SIZE=${DEFAULT_MESH_SIZE})

add_executable(exercise-01 src/exercise-01.cpp src/args.cpp src/Current.cpp)
deal_ii_setup_target(exercise-01)
