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
add_subdirectory(external/repair_gn)

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
        ${CMAKE_HOME_DIRECTORY}/external/repair_gn
        ${STXXL_INCLUDE_DIRS}
        )


add_executable(estimation_size src/estimation_size.cpp)
target_link_libraries(estimation_size sdsl divsufsort divsufsort64)

add_executable(rct_index src/rct_main.cpp)
target_link_libraries(rct_index sdsl divsufsort divsufsort64)

add_executable(rct_index_gn src/rct_main_gn.cpp)
target_link_libraries(rct_index_gn sdsl divsufsort divsufsort64)

add_executable(rct_index_gn_redundant src/rct_main_gn_redundant.cpp)
target_link_libraries(rct_index_gn_redundant sdsl divsufsort divsufsort64)

add_executable(rct_index_gn_separately src/rct_main_gn_separately.cpp)
target_link_libraries(rct_index_gn_separately sdsl divsufsort divsufsort64)

add_executable(rct_index_gn_max_wr_min src/rct_main_gn_max_wr_min.cpp)
target_link_libraries(rct_index_gn_max_wr_min sdsl divsufsort divsufsort64)

add_executable(rct_index_dac src/rct_dac_main.cpp)
target_link_libraries(rct_index_dac sdsl divsufsort divsufsort64)

add_executable(alternative_code src/alternative_code.cpp)
target_link_libraries(alternative_code sdsl divsufsort divsufsort64)

add_executable(rlz_test src/rlz_test.cpp)
target_link_libraries(rlz_test sdsl divsufsort divsufsort64)

add_executable(rlz src/rlz.cpp)
target_link_libraries(rlz sdsl divsufsort divsufsort64)

add_executable(rlz_file src/rlz_file.cpp)
target_link_libraries(rlz_file sdsl divsufsort divsufsort64)

add_executable(stats_rct src/stats_rct.cpp)
target_link_libraries(stats_rct sdsl divsufsort divsufsort64)

add_executable(runs_bitvector_test src/runs_bitvector_test.cpp)
target_link_libraries(runs_bitvector_test sdsl divsufsort divsufsort64)

add_executable(cut_dataset src/cut_dataset.cpp)
target_link_libraries(cut_dataset sdsl divsufsort divsufsort64)

add_executable(rct_test src/rct_test.cpp)
target_link_libraries(rct_test sdsl divsufsort divsufsort64)

add_executable(rct_gn_test src/rct_gn_test.cpp)
target_link_libraries(rct_gn_test sdsl divsufsort divsufsort64)

add_executable(rct_test_repair src/rct_test_repair.cpp)
target_link_libraries(rct_test_repair sdsl divsufsort divsufsort64)

add_executable(rct_times_multiple src/rct_times_multiple.cpp)
target_link_libraries(rct_times_multiple sdsl divsufsort divsufsort64)

add_executable(rct_times_multiple_sparse src/rct_times_multiple_sparse.cpp)
target_link_libraries(rct_times_multiple_sparse sdsl divsufsort divsufsort64)

add_executable(rct_times src/rct_times.cpp)
target_link_libraries(rct_times sdsl divsufsort divsufsort64)

add_executable(rct_times_gn src/rct_times_gn.cpp)
target_link_libraries(rct_times_gn sdsl divsufsort divsufsort64)

add_executable(rct_times_gn_max_wr_min src/rct_times_gn_max_wr_min.cpp)
target_link_libraries(rct_times_gn_max_wr_min sdsl divsufsort divsufsort64)

add_executable(rct_times_gn_redundant src/rct_times_gn_redundant.cpp)
target_link_libraries(rct_times_gn_redundant sdsl divsufsort divsufsort64)

add_executable(rct_times_gn_separately src/rct_times_gn_separately.cpp)
target_link_libraries(rct_times_gn_separately sdsl divsufsort divsufsort64)

add_executable(rct_fix_disap src/rct_fix_disap.cpp)
target_link_libraries(rct_fix_disap sdsl divsufsort divsufsort64)

add_executable(rct_index_rtree src/rct_main_v2.cpp)
target_link_libraries(rct_index_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(v1_to_v2 src/v1_to_v2.cpp)
target_link_libraries(v1_to_v2 sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(v1_to_v2_repair src/v1_to_v2_repair.cpp)
target_link_libraries(v1_to_v2_repair sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(v1_to_v2_multiple src/v1_to_v2_multiple.cpp)
target_link_libraries(v1_to_v2_multiple sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_test_repair_rtree src/rct_test_repair_rtree.cpp)
target_link_libraries(rct_test_repair_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_test_multiple_rtree src/rct_test_multiple_rtree.cpp)
target_link_libraries(rct_test_multiple_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_test_sparse_multiple_rtree src/rct_test_multiple_sparse_rtree.cpp)
target_link_libraries(rct_test_sparse_multiple_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_test_sparse_multiple_rtree_no_min_max src/rct_test_multiple_sparse_rtree_no_min_max.cpp)
target_link_libraries(rct_test_sparse_multiple_rtree_no_min_max sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_test_multiple_rtree_no_min_max src/rct_test_multiple_rtree_no_min_max.cpp)
target_link_libraries(rct_test_multiple_rtree_no_min_max sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_test_multiple_rtree_lite src/rct_test_multiple_rtree_lite.cpp)
target_link_libraries(rct_test_multiple_rtree_lite sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(to_no_min_max src/to_no_min_max.cpp)
target_link_libraries(to_no_min_max sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(to_lite src/to_lite.cpp)
target_link_libraries(to_lite sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(to_lite_sparse src/to_lite_sparse.cpp)
target_link_libraries(to_lite_sparse sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(to_no_min_max_sparse src/to_no_min_max_sparse.cpp)
target_link_libraries(to_no_min_max_sparse sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_test_rtree src/rct_test_rtree.cpp)
target_link_libraries(rct_test_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rtree_test src/rtree_test.cpp)
target_link_libraries(rtree_test sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_times_rtree src/rct_times_rtree.cpp)
target_link_libraries(rct_times_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_times_repair_rtree src/rct_times_repair_rtree.cpp)
target_link_libraries(rct_times_repair_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_times_repair src/rct_times_repair.cpp)
target_link_libraries(rct_times_repair sdsl divsufsort divsufsort64)

add_executable(rct_times_multiple_rtree src/rct_times_multiple_rtree.cpp)
target_link_libraries(rct_times_multiple_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_times_multiple_sparse_rtree src/rct_times_multiple_sparse_rtree.cpp)
target_link_libraries(rct_times_multiple_sparse_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_times_multiple_rtree_no_min_max src/rct_times_multiple_rtree_no_min_max.cpp)
target_link_libraries(rct_times_multiple_rtree_no_min_max sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_times_multiple_sparse_rtree_no_min_max src/rct_times_multiple_sparse_rtree_no_min_max.cpp)
target_link_libraries(rct_times_multiple_sparse_rtree_no_min_max sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_times_multiple_rtree_lite src/rct_times_multiple_rtree_lite.cpp)
target_link_libraries(rct_times_multiple_rtree_lite sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_times_multiple_sparse_rtree_lite src/rct_times_multiple_sparse_rtree_lite.cpp)
target_link_libraries(rct_times_multiple_sparse_rtree_lite sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(rct_repair src/rct_main_repair.cpp)
target_link_libraries(rct_repair libRepairGN sdsl divsufsort divsufsort64 spatialindex STRRtree)

#add_executable(rct_repair_v2 src/rct_main_repair_v2.cpp)
#target_link_libraries(rct_repair_v2 libRepairGN sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(lz_parser_repair_test src/lz_parser_repair_test.cpp)
target_link_libraries(lz_parser_repair_test  libRepairGN sdsl divsufsort divsufsort64)

add_executable(sa src/sa.cpp)
target_link_libraries(sa sdsl divsufsort divsufsort64)

add_executable(rct_multiple src/rct_main_multiple.cpp)
target_link_libraries(rct_multiple sdsl divsufsort divsufsort64)

add_executable(rct_multiple_rev src/rct_main_multiple_rev.cpp)
target_link_libraries(rct_multiple_rev sdsl divsufsort divsufsort64)

add_executable(rct_multiple_test src/rct_test_multiple.cpp)
target_link_libraries(rct_multiple_test sdsl divsufsort divsufsort64)

add_executable(rct_multiple_sparse_test src/rct_test_multiple_sparse.cpp)
target_link_libraries(rct_multiple_sparse_test sdsl divsufsort divsufsort64)

add_executable(update_period src/update_period.cpp)
target_link_libraries(update_period sdsl divsufsort divsufsort64)

add_executable(update_period_repair src/update_period_repair.cpp)
target_link_libraries(update_period_repair sdsl divsufsort divsufsort64)

add_executable(update_period_multiple src/update_period_multiple.cpp)
target_link_libraries(update_period_multiple sdsl divsufsort divsufsort64)

add_executable(size_index src/size_index.cpp)
target_link_libraries(size_index sdsl divsufsort divsufsort64)

add_executable(size_index_rtree src/size_index_rtree.cpp)
target_link_libraries(size_index_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(size_multiple_index src/size_multiple_index.cpp)
target_link_libraries(size_multiple_index sdsl divsufsort divsufsort64)

add_executable(size_multiple_index_rtree src/size_multiple_index_rtree.cpp)
target_link_libraries(size_multiple_index_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(size_multiple_sparse_index_rtree src/size_multiple_sparse_index_rtree.cpp)
target_link_libraries(size_multiple_sparse_index_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(size_repair_index src/size_repair_index.cpp)
target_link_libraries(size_repair_index sdsl divsufsort divsufsort64)

add_executable(size_repair_index_rtree src/size_repair_index_rtree.cpp)
target_link_libraries(size_repair_index_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(sparse_multiple src/sparse_multiple.cpp)
target_link_libraries(sparse_multiple sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(sparse_multiple_rtree src/sparse_multiple_rtree.cpp)
target_link_libraries(sparse_multiple_rtree sdsl divsufsort divsufsort64 spatialindex STRRtree)

add_executable(proba test/test.cpp)
target_link_libraries(proba sdsl divsufsort divsufsort64)

add_executable(loading_reference src/loading_reference.cpp)
target_link_libraries(loading_reference sdsl divsufsort divsufsort64)

add_executable(baseline_times src/baseline_times.cpp)
target_link_libraries(baseline_times sdsl divsufsort divsufsort64)

add_executable(space_estimation src/space_estimation.cpp)
target_link_libraries(space_estimation sdsl divsufsort divsufsort64 spatialindex STRRtree)
