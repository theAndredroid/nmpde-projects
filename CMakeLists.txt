cmake_minimum_required(VERSION 3.12.0)
project(project_title LANGUAGES CXX C)

include(cmake-common.cmake)

# add_executable(exercise-01 src/exercise-01.cpp src/Current.cpp)
add_executable(exercise-01 src/Current.cpp)
deal_ii_setup_target(exercise-01)
