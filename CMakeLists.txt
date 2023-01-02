cmake_minimum_required (VERSION 2.8.12)
project (alignment-writer)

set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Determine build type and set flags accordingly
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) ## Build Release by default
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpedantic -g")

## Set C++17 support depending on cmake version
if (${CMAKE_MAJOR_VERSION} GREATER 2 AND ${CMAKE_MINOR_VERSION} GREATER 0)
  set (CMAKE_CXX_STANDARD 17)
  set (CMAKE_CXX_STANDARD_REQUIRED ON)
else()
  add_compile_options(-std=c++17)
endif()

## alignment-writer library
add_library(libalignmentwriter
  ${CMAKE_CURRENT_SOURCE_DIR}/src/unpack.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/pack.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/merge.cpp)
set_target_properties(libalignmentwriter PROPERTIES OUTPUT_NAME alignment-writer)

## alignment-writer executable
add_executable(alignment-writer ${CMAKE_CURRENT_SOURCE_DIR}/src/alignment_writer.cpp)
target_link_libraries(alignment-writer libalignmentwriter)

## Own header files
set(CMAKE_ALIGNMENT_WRITER_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_ALIGNMENT_WRITER_HEADERS})

## Check supported compression types
find_package(BZip2)
if (BZIP2_FOUND)
  include_directories(${BZIP2_INCLUDE_DIRS})
  target_link_libraries(alignment-writer ${BZIP2_LIBRARIES})
endif()
find_package(LibLZMA)
if (LIBLZMA_FOUND)
  include_directories(${LIBLZMA_INCLUDE_DIRS})
  target_link_libraries(alignment-writer ${LIBLZMA_LIBRARIES})
endif()
find_package(ZLIB)
if (ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(alignment-writer ${ZLIB_LIBRARIES})
endif()

## Check dependencies and download them if not found
## BitMagic - supplied with the project
set(CMAKE_BITMAGIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/BitMagic-7.12.3/src)
include_directories(${CMAKE_BITMAGIC_HEADERS})

## bxzstr
if (DEFINED CMAKE_BXZSTR_HEADERS)
  message(STATUS "bxzstr headers provided in: ${CMAKE_BXZSTR_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-bxzstr.txt.in ${CMAKE_BINARY_DIR}/external/bxzstr-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "CMake step for bxzstr failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/bxzstr-download )
  if(result)
    message(FATAL_ERROR "Build step for bxzstr failed: ${result}")
  endif()
  set(CMAKE_BXZSTR_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/bxzstr/include)
endif()
include_directories(${CMAKE_BXZSTR_HEADERS})

## cxxargs
if (DEFINED CMAKE_CXXARGS_HEADERS)
  message(STATUS "cxxargs headers provided in: ${CMAKE_CXXARGS_HEADERS}")
else()
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/CMakeLists-cxxargs.txt.in ${CMAKE_BINARY_DIR}/external/cxxargs-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/cxxargs-download )
  if(result)
    message(FATAL_ERROR "CMake step for cxxargs failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/external/cxxargs-download )
  if(result)
    message(FATAL_ERROR "Build step for cxxargs failed: ${result}")
  endif()
  set(CMAKE_CXXARGS_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/external/cxxargs/include)
endif()
include_directories(${CMAKE_CXXARGS_HEADERS})

## Include project headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

## Get version number from git tags
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_VARIABLE ALIGNMENT_WRITER_BUILD_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message( STATUS "alignment-writer version: ${ALIGNMENT_WRITER_BUILD_VERSION}")
else()
  message(STATUS "GIT not found, version numbering disabled")
  set(ALIGNMENT_WRITER_BUILD_VERSION "unknown")
endif()
string(TIMESTAMP ALIGNMENT_WRITER_BUILD_TIMESTAMP)

## Generate a version.h file containing build version and timestamp
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/version.h @ONLY)
