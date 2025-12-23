
cmake_minimum_required(VERSION 3.22)
project(radt LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable verbose output
set(CMAKE_VERBOSE_MAKEFILE ON)

# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# HYPRE installation path
set(HYPRE_DIR "/home/mike/proj/hypre/src/hypre")
include_directories(${HYPRE_DIR}/include)
link_directories(${HYPRE_DIR}/lib)

find_library(HYPRE_LIB NAMES HYPRE PATHS ${HYPRE_DIR}/lib NO_DEFAULT_PATH)

if(NOT HYPRE_LIB)
    message(FATAL_ERROR "Could not find HYPRE library in ${HYPRE_DIR}/lib")
endif()
message(STATUS "Found HYPRE library: ${HYPRE_LIB}")

# Include your own headers
include_directories(include)

# Source files
set(SOURCES
    src/geometry.cpp
    src/Material.cpp
    src/physics.cpp
    src/solver.cpp
    src/setup.cpp
    src/utilities.cpp
)

# Create executable
add_library(radt_lib ${SOURCES})
target_include_directories(radt_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(radt src/main.cpp)

# Link libraries
target_link_libraries(radt ${MPI_CXX_LIBRARIES} radt_lib ${HYPRE_LIB})
add_subdirectory(tests)
