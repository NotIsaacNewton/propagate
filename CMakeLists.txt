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

add_executable(propagate propagate.cpp filetools.cpp fftw_complex_tools.cpp console_tools.cpp interpolate_1d.cpp)
target_link_libraries(propagate ${FFT_LIBS})

add_executable(wavepacket wavepackets.cpp filetools.cpp fftw_complex_tools.cpp console_tools.cpp interpolate_1d.cpp)
target_link_libraries(wavepacket ${FFT_LIBS})

add_executable(potential potentials.cpp filetools.cpp fftw_complex_tools.cpp console_tools.cpp interpolate_1d.cpp)
target_link_libraries(potential ${FFT_LIBS})

add_executable(interpolate interpolate_1d.cpp filetools.cpp console_tools.cpp)
target_link_libraries(interpolate)
