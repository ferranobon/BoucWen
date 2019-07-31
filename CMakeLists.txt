# Name of the project
cmake_minimum_required(VERSION 3.12.0 FATAL_ERROR)
project (BoucWen)

# Build type
if( NOT (CMAKE_BUILD_TYPE MATCHES Debug ) AND
    NOT (CMAKE_BUILD_TYPE MATCHES Release ))
  message( FATAL_ERROR
    "CMAKE_BUILD_TYPE must be one of: Debug or Release.")
endif()

# Set variables for debug
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(DEBUG TRUE)
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(RELEASE TRUE)
endif()

# Make will display verbose output
set(CMAKE_VERBOSE_MAKEFILE ON)

# Compile options
option(PROFILING      "Enable profiling using gprof"            OFF)

if(DEBUG)
  # -D_POSIX_SOURCE is added in order to make gcc be POSIX compliant
  if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "-m64 -std=c++11 -pedantic -Wall -W -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wbad-function-cast -Wundef -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Wstrict-overflow=5 -fshort-enums -fno-common -g3")
  endif()
  if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -stdlib=libc++ -pedantic -Wall -Wextra -Wmissing-prototypes -Wstrict-prototypes -Wbad-function-cast -Wundef -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wnested-externs -Wstrict-overflow=5 -fshort-enums -fno-common -g")
  endif()
  if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -pedantic -Wall -Wmissing-prototypes -Wstrict-prototypes -Wshadow -Wpointer-arith -Wwrite-strings -debug all" )
  endif()
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS_DEBUG "-Minform=inform")
  endif()
  # Add profiling flags.
  if(PROFILING)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg") 
 endif() 
 list(APPEND COMPILE_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})
endif()

if(RELEASE)
  set(CMAKE_CXX_FLAGS_RELEASE "-march=native -O2")
  # Add profiling flags.
  if(PROFILING)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -pg")
  endif()
  list(APPEND COMPILE_FLAGS ${CMAKE_C_FLAGS_RELEASE})
endif()

# Find Boost
find_package(Boost 1.56 REQUIRED COMPONENTS program_options)

add_executable(BoucWen main.cpp)
target_link_libraries(BoucWen Boost::program_options)

