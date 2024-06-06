cmake_minimum_required (VERSION 3.11)
project (alignment-writer)
include (FetchContent)

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
  ${CMAKE_CURRENT_SOURCE_DIR}/src/pack.cpp)
set_target_properties(libalignmentwriter PROPERTIES OUTPUT_NAME alignment-writer)

## alignment-writer executable
add_executable(alignment-writer ${CMAKE_CURRENT_SOURCE_DIR}/src/alignment_writer.cpp)
target_link_libraries(alignment-writer libalignmentwriter)

## Check dependencies
find_package(OpenMP)
if (OPENMP_FOUND)
  set(ALIGNMENTWRITER_OPENMP_SUPPORT 1)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_GLIBCXX_PARALLEL")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -D_GLIBCXX_PARALLEL")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -D_GLIBCXX_PARALLEL")
else()
  set(ALIGNMENTWRITER_OPENMP_SUPPORT 0)
endif()

## Configure OpenMP if it is supported on the system.
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/alignment-writer_openmp_config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/include/alignment-writer_openmp_config.hpp @ONLY)

## Own header files
set(CMAKE_ALIGNMENT_WRITER_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_ALIGNMENT_WRITER_HEADERS})

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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/include/version.h @ONLY)

## Check supported compression types
#### zlib
if ((DEFINED ZLIB_LIBRARY AND DEFINED ZLIB_INCLUDE_DIR) AND (NOT DEFINED ZLIB_FOUND))
  message(STATUS "zlib library provided in: " ${ZLIB_LIBRARY})
  message(STATUS "zlib headers provided in: " ${ZLIB_INCLUDE_DIR})
  include_directories(${ZLIB_INCLUDE_DIR})
  target_link_libraries(alignment-writer ${ZLIB_LIBRARY})
  set(ALIGNMENT_WRITER_HAVE_ZLIB 1)
else()
  find_package(ZLIB)
  if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIR})
    target_link_libraries(alignment-writer ${ZLIB_LIBRARY})
    set(ALIGNMENT_WRITER_HAVE_ZLIB 1)
  else()
    set(ALIGNMENT_WRITER_HAVE_ZLIB 0)
  endif()
endif()

#### bzip2
if (DEFINED BZIP2_LIBRARIES AND DEFINED BZIP2_INCLUDE_DIR AND (NOT DEFINED BZIP2_FOUND))
  message(STATUS "bzip2 library provided in: " ${BZIP2_LIBRARIES})
  message(STATUS "bzip2 headers provided in: " ${BZIP2_INCLUDE_DIR})
  include_directories(${BZIP2_INCLUDE_DIR})
  target_link_libraries(alignment-writer ${BZIP2_LIBRARIES})
  set(ALIGNMENT_WRITER_HAVE_BZIP2 1)
else()
  find_package(BZip2)
  if (BZIP2_FOUND)
    include_directories(${BZIP2_INCLUDE_DIR})
    target_link_libraries(alignment-writer ${BZIP2_LIBRARIES})
    set(ALIGNMENT_WRITER_HAVE_BZIP2 1)
  else()
    set(ALIGNMENT_WRITER_HAVE_BZIP2 0)
  endif()
endif()

#### lzma
if (DEFINED LIBLZMA_LIBRARY AND DEFINED LIBLZMA_INCLUDE_DIR AND (NOT DEFINED LIBLZMA_FOUND))
  message(STATUS "liblzma library provided in: " ${LIBLZMA_LIBRARY})
  message(STATUS "liblzma headers provided in: " ${LIBLZMA_INCLUDE_DIR})
  include_directories(${LIBLZMA_INCLUDE_DIR})
  target_link_libraries(alignment-writer ${LIBLZMA_LIBRARY})
  set(ALIGNMENT_WRITER_HAVE_LIBLZMA 1)
else()
  find_package(LibLZMA)
  if (LIBLZMA_FOUND)
    include_directories(${LIBLZMA_INCLUDE_DIR})
    target_link_libraries(alignment-writer ${LIBLZMA_LIBRARY})
    set(ALIGNMENT_WRITER_HAVE_LIBLZMA 1)
  else()
    set(ALIGNMENT_WRITER_HAVE_LIBLZMA 0)
  endif()
endif()

## Check dependencies and download them if not found
### BitMagic - supplied with the project
set(CMAKE_BITMAGIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/BitMagic-7.12.3/src)
include_directories(${CMAKE_BITMAGIC_HEADERS})

### kseq++
set(CMAKE_KSEQPP_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/kseqpp-1.1.2/include)
include_directories(${CMAKE_KSEQPP_HEADERS})

## bxzstr
if (DEFINED CMAKE_BXZSTR_HEADERS)
  message(STATUS "bxzstr headers provided in: ${CMAKE_BXZSTR_HEADERS}")
else()
  FetchContent_Declare(bxzstr
    GIT_REPOSITORY    https://github.com/tmaklin/bxzstr.git
    GIT_TAG           v1.1.0
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/bxzstr"
    BUILD_IN_SOURCE   0
    CMAKE_ARGS        -D ZSTD_FOUND=0
                      -D ZLIB_FOUND=${ALIGNMENT_WRITER_HAVE_ZLIB}
		      -D BZIP2_FOUND=${ALIGNMENT_WRITER_HAVE_BZIP2}
		      -D LIBLZMA_FOUND=${ALIGNMENT_WRITER_HAVE_LIBLZMA}
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(bxzstr)
  set(CMAKE_BXZSTR_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/bxzstr/include)
endif()
include_directories(${CMAKE_BXZSTR_HEADERS})

## cxxargs
if (DEFINED CMAKE_CXXARGS_HEADERS)
  message(STATUS "cxxargs headers provided in: ${CMAKE_CXXARGS_HEADERS}")
else()
  FetchContent_Declare(cxxargs
    GIT_REPOSITORY    https://github.com/tmaklin/cxxargs.git
    GIT_TAG           v1.1.4
    PREFIX            "external"
    SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/external/cxxargs"
    BUILD_IN_SOURCE   0
    BUILD_COMMAND     ""
    CONFIGURE_COMMAND ""
    INSTALL_COMMAND   ""
    )
  FetchContent_MakeAvailable(cxxargs)
  set(CMAKE_CXXARGS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/external/cxxargs/include)
endif()
include_directories("${CMAKE_CXXARGS_HEADERS}")

# Link libraries
if (OPENMP_FOUND)
  target_link_libraries(libalignmentwriter OpenMP::OpenMP_CXX)
endif()
