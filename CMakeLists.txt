cmake_minimum_required (VERSION 3.7)
project (cglass
    VERSION 0.1.0
    DESCRIPTION "A C++ project using CMake"
    LANGUAGES CXX
)

## Require C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/.CMake_Modules")

# Create a version number for the software from the git sha1 of the current commit
# Can be referenced as "@GIT_SHA1@" in a header file from a compiler definition
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
#add_definitions(-DSC_VERSION="${GIT_SHA1}")

if(DEBUG OR TRACE)
  set(CMAKE_BUILD_TYPE "Debug")
else()
  set(CMAKE_BUILD_TYPE "Release")
endif()

if (OMP AND APPLE)
  set(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
  set(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")
  set(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
  set(OPENMP_INCLUDES "/usr/local/opt/llvm/include")
endif()

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

## Add external dependencies here, e.g.:
add_subdirectory(extern/KMC)

## The directory for the library
add_subdirectory(src)

## Only build documentation if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    ## Docs only available if this is the main project
    find_package(Doxygen QUIET)
    if(Doxygen_FOUND)
        add_subdirectory(docs)
    else()
        message(STATUS "Doxygen not found, not building docs")
    endif()
    ## The example executable is only available for main project
    #add_subdirectory(examples)
endif()

## Allow testing for main and subprojects if -DTESTS=TRUE
if(TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
