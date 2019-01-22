cmake_minimum_required(VERSION 2.8.7)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)

project(RCT)
set(PROJECT_VENDOR "Adrián Gómez Brandón")
set(CMAKE_PROJECT_NAME "RCT")
set(PROJECT_CONTACT "adrian.gbrandon@udc.es")
set(CMAKE_INSTALL_PREFIX /home/adrian/software)

set(CMAKE_BUILD_TYPE "Release")

# C++11 compiler Check
if(NOT CMAKE_CXX_COMPILER_VERSION) # work around for cmake versions smaller than 2.8.10
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE CMAKE_CXX_COMPILER_VERSION)
endif()
if(CMAKE_CXX_COMPILER MATCHES ".*clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    set(CMAKE_COMPILER_IS_CLANGXX 1)
endif()
if( (CMAKE_COMPILER_IS_GNUCXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 4.7) OR
(CMAKE_COMPILER_IS_CLANGXX AND ${CMAKE_CXX_COMPILER_VERSION} VERSION_LESS 3.2))
    message(FATAL_ERROR "Your C++ compiler does not support C++11. Please install g++ 4.7 (or greater) or clang 3.2 (or greater)")
else()
    message(STATUS "Compiler is recent enough to support C++11.")
endif()
if( CMAKE_COMPILER_IS_GNUCXX )
    append_cxx_compiler_flags("-std=c++11 -Wall -Wextra " "GCC" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-msse4.2 -O3 -ffast-math -funroll-loops -fno-omit-frame-pointer -g" "GCC" CMAKE_CXX_FLAGS_RELEASE)
else()
    append_cxx_compiler_flags("-std=c++11" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-stdlib=libc++" "CLANG" CMAKE_CXX_FLAGS)
    append_cxx_compiler_flags("-msse4.2 -O3  -ffast-math -funroll-loops -DNDEBUG" "CLANG" CMAKE_CXX_FLAGS_RELEASE)
endif()
include(CheckSSE)
FindSSE ()
if( SSE4_2_FOUND )
    if( CMAKE_COMPILER_IS_GNUCXX )
        append_cxx_compiler_flags("-msse4.2" "GCC" CMAKE_CXX_FLAGS)
    else()
        append_cxx_compiler_flags("-msse4.2" "CLANG" CMAKE_CXX_FLAGS)
    endif()
    message(STATUS "CPU does support SSE4.2.")
else()
    message(STATUS "CPU does NOT support SSE4.2")
endif()

add_subdirectory(external/sdsl-lite)
add_subdirectory(external/STRRtree)

include_directories("${PROJECT_SOURCE_DIR}/external/googletest/include"
        ${CMAKE_HOME_DIRECTORY}/include
        ${CMAKE_BINARY_DIR}/external/sdsl-lite/include
        ${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort/include/
        ${CMAKE_HOME_DIRECTORY}/external/STRRtree/include
        ${CMAKE_HOME_DIRECTORY}/external/STRRtree/rdt_wt
        ${CMAKE_HOME_DIRECTORY}/external/STRRtree/util
        ${CMAKE_HOME_DIRECTORY}/external/STRRtree/external/libspatialindex/include/
        ${CMAKE_HOME_DIRECTORY}/external/STRRtree/external/libspatialindex/src/
        ${CMAKE_HOME_DIRECTORY}/external/googletest/include
        ${CMAKE_HOME_DIRECTORY}/external/sdsl-lite/external/FastPFor/headers/
        ${STXXL_INCLUDE_DIRS}
        )


add_executable(rlz_file src/rlz_file.cpp)
target_link_libraries(rlz_file sdsl divsufsort divsufsort64)

add_executable(estimation_size src/estimation_size.cpp)
target_link_libraries(estimation_size sdsl divsufsort divsufsort64)

add_executable(rct_index src/rct_main.cpp include/runs_vector.hpp)
target_link_libraries(rct_index sdsl divsufsort divsufsort64)

add_executable(rct_index_dac src/rct_dac_main.cpp)
target_link_libraries(rct_index_dac sdsl divsufsort divsufsort64)

add_executable(alternative_code src/alternative_code.cpp)
target_link_libraries(alternative_code sdsl divsufsort divsufsort64)

add_executable(rlz_test src/rlz_test.cpp)
target_link_libraries(rlz_test sdsl divsufsort divsufsort64)

add_executable(stats_rct src/stats_rct.cpp)
target_link_libraries(stats_rct sdsl divsufsort divsufsort64)

add_executable(runs_bitvector_test src/runs_bitvector_test.cpp)
target_link_libraries(runs_bitvector_test sdsl divsufsort divsufsort64)

add_executable(cut_dataset src/cut_dataset.cpp)
target_link_libraries(cut_dataset sdsl divsufsort divsufsort64)

add_executable(rct_test src/rct_test.cpp)
target_link_libraries(rct_test sdsl divsufsort divsufsort64)

add_executable(rct_times src/rct_times.cpp)
target_link_libraries(rct_times sdsl divsufsort divsufsort64)

add_executable(rct_fix_disap src/rct_fix_disap.cpp)
target_link_libraries(rct_fix_disap sdsl divsufsort divsufsort64)

add_executable(rct_index_rtree src/rct_main_v2.cpp)
target_link_libraries(rct_index_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(v1_to_v2 src/v1_to_v2.cpp)
target_link_libraries(v1_to_v2 sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(proba src/proba.cpp)

add_executable(rct_test_v2 src/rct_test_v2.cpp)
target_link_libraries(rct_test_v2 sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rtree_test src/rtree_test.cpp)
target_link_libraries(rtree_test sdsl divsufsort divsufsort64 spatialindex STRRtree)