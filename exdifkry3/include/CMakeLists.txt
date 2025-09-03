cmake_minimum_required(VERSION 3.10)

project(RadiativeTransfer2D)

 

# Set C++ standard

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

 

# Enable verbose output

set(CMAKE_VERBOSE_MAKEFILE ON)

 

# Find MPI

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

 

# HYPRE installation path (update this if needed)

set(HYPRE_DIR "/usr/local/hypre")

include_directories(${HYPRE_DIR}/include)

link_directories(${HYPRE_DIR}/lib)

 

# Source files

file(GLOB SRC_FILES

    src/*.cpp

    src/*.h

)

 

# Executable

add_executable(radiative_transfer ${SRC_FILES})

 

# Link libraries

target_link_libraries(radiative_transfer

    ${MPI_CXX_LIBRARIES}

    HYPRE

)

 

# Optional: show compile flags

message(STATUS "MPI Include Path: ${MPI_INCLUDE_PATH}")

message(STATUS "HYPRE Include Path: ${HYPRE_DIR}/include")