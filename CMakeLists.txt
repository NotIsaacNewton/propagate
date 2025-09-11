cmake_minimum_required(VERSION 3.29)
project(propagate)

set(CMAKE_CXX_STANDARD 23)

#option(BUILD_INTEL "Build for Beocat x86_64" OFF)
#if(BUILD_INTEL)
#    set(CMAKE_OSX_ARCHITECTURES "x86_64")
#    add_compile_options(-march=x86-64 -mtune=sandybridge)
#endif()

set(ARMPL_DIR /opt/arm/armpl_25.04_flang-new_clang_19)

include_directories(${ARMPL_DIR}/include)
link_directories(${ARMPL_DIR}/lib)

add_executable(propagate propagate.cpp)
target_link_libraries(propagate armpl)

add_executable(wavepacket wavepackets.cpp)
target_link_libraries(wavepacket armpl)

add_executable(potential potentials.cpp)
target_link_libraries(potential armpl)
