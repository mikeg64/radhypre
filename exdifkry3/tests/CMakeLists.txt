
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
set(GOOGLETEST_DIR "/home/mike/proj/googletest/build")
include_directories(${HYPRE_DIR}/include
                    ${GOOGLETEST_DIR}/include)
link_directories(${HYPRE_DIR}/lib
                 ${GOOGLETEST_DIR}/lib)
                 

# Include your own headers
include_directories(include)

# Source files
set(SOURCES
    ../src/main.cpp
    ../src/geometry.cpp
    ../src/Material.cpp
    ../src/physics.cpp
    ../src/solver.cpp
    ../src/setup.cpp
    ../src/utilities.cpp
    hyptest1.cpp
    hyptest1_unittest.cpp
    hyptest
)

# Create executable
add_executable(hyptest ${SOURCES})

# Link libraries
target_link_libraries(hyptest ${MPI_CXX_LIBRARIES} HYPRE)
