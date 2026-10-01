cmake_minimum_required(VERSION 3.29)
project(propagate)

set(CMAKE_CXX_STANDARD 23)

# flag to build for intel chips
option(BUILD_INTEL "Build for x86_64 using real FFTW3" OFF)

# configure OpenMP
if(NOT BUILD_INTEL)
    set(OpenMP_ROOT "/opt/homebrew/opt/libomp")
endif()
find_package(OpenMP REQUIRED)

# configure fftw3
if(BUILD_INTEL)
    set(FFT_LIBS fftw3 m)
else()
    set(ARMPL_DIR /opt/arm/armpl_25.04_flang-new_clang_19)
    include_directories(${ARMPL_DIR}/include)
    link_directories(${ARMPL_DIR}/lib)
    set(FFT_LIBS armpl)
endif()

# find, add, and target library of shared tools
file(GLOB_RECURSE lib CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp")
add_library(LIB STATIC ${lib})
target_include_directories(LIB PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib ${OpenMP_ROOT}/include)
target_link_libraries(LIB PUBLIC OpenMP::OpenMP_CXX)

# find directories with main programs and source codes
set(apps ${CMAKE_CURRENT_SOURCE_DIR}/apps)
set(propagate ${apps}/propagate)
set(potential ${apps}/potentials)
set(wavepacket ${apps}/wavepackets)
file(GLOB_RECURSE PROPAGATE CONFIGURE_DEPENDS "${propagate}/*.cpp")
file(GLOB_RECURSE POTENTIAL CONFIGURE_DEPENDS "${potential}/*.cpp")
file(GLOB_RECURSE WAVEPACKET CONFIGURE_DEPENDS "${wavepacket}/*.cpp")

# compile with optimizations
target_compile_options(LIB PUBLIC -O3 -march=native)

# build propagate
add_executable(propagate ${PROPAGATE})
target_link_libraries(propagate LIB ${FFT_LIBS} OpenMP::OpenMP_CXX)

# build wavepacket
add_executable(wavepacket ${WAVEPACKET})
target_link_libraries(wavepacket LIB ${FFT_LIBS} OpenMP::OpenMP_CXX)

# build potential
add_executable(potential ${POTENTIAL})
target_link_libraries(potential LIB ${FFT_LIBS} OpenMP::OpenMP_CXX)
