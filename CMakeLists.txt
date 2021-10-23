cmake_minimum_required(VERSION 3.2)
# project and version
project(bam_stats LANGUAGES CXX VERSION 0.1.0)

# Standard GNU install dirs
# include(GNUInstallDirs)

# PROJECT_BINARY_DIR default is the build dir
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# CXX configure
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -pedantic")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set static library
# set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build all libraries as static")

# htslib library and include files in vendor
set(HTS_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib/htslib)
set(HTS_LIB ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib/libhts.a)

if (NOT EXISTS ${HTS_INCLUDE})
    message(FATAL_ERROR "[FATAL_ERROR]: Can not find htslib directory in \
    ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib, Please check vendor directory \
    for existence of htslib.")
endif (NOT EXISTS ${HTS_INCLUDE})

if (EXISTS ${HTS_LIB})
    message(STATUS "[STATUS]: Find htslib static library libhts.a")
else (EXISTS ${HTS_LIB})
    message(STATUS "[STATUS]: htslib static library libhts.a not exist in ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib, build one ...")
    execute_process(COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib)
endif (EXISTS ${HTS_LIB})

add_executable(bam_stats
    src/bam_stats.cpp
    src/from_cigar.cpp)

target_include_directories(bam_stats PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include>)

target_link_libraries(bam_stats ${CMAKE_CURRENT_SOURCE_DIR}/vendor/htslib/libhts.a z m bz2 lzma pthread curl)

