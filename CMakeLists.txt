cmake_minimum_required(VERSION 3.29)
project(propagate)

set(CMAKE_CXX_STANDARD 23)

set(ARMPL_DIR /opt/arm/armpl_25.04_flang-new_clang_19)

include_directories(${ARMPL_DIR}/include)
link_directories(${ARMPL_DIR}/lib)

add_executable(propagate main.cpp)

target_link_libraries(propagate armpl)
