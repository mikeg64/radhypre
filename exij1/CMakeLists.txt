cmake_minimum_required(VERSION 3.10)
project(hypre_example CXX)

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




add_executable(hypre_example main.cpp)
target_link_libraries(hypre_example MPI::MPI_CXX ${HYPRE_LIB})