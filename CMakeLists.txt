cmake_minimum_required(VERSION 3.29)
project(propagate)

set(CMAKE_CXX_STANDARD 23)

option(BUILD_INTEL "Build for x86_64 (homelab/Beocat) using real FFTW3" OFF)

# option to build for intel chips. for HPCs, it may be necessary to do a static link.
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
target_include_directories(LIB PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# find directories with main programs and source codes
set(apps ${CMAKE_CURRENT_SOURCE_DIR}/apps)
set(propagate ${apps}/propagate)
set(potential ${apps}/potentials)
set(wavepacket ${apps}/wavepackets)
file(GLOB_RECURSE PROPAGATE CONFIGURE_DEPENDS "${propagate}/*.cpp")
file(GLOB_RECURSE POTENTIAL CONFIGURE_DEPENDS "${potential}/*.cpp")
file(GLOB_RECURSE WAVEPACKET CONFIGURE_DEPENDS "${wavepacket}/*.cpp")

# build propagate
add_executable(propagate ${PROPAGATE})
target_link_libraries(propagate LIB ${FFT_LIBS})

# build wavepacket
add_executable(wavepacket ${WAVEPACKET})
target_link_libraries(wavepacket LIB ${FFT_LIBS})

# build potential
add_executable(potential ${POTENTIAL})
target_link_libraries(potential LIB ${FFT_LIBS})
