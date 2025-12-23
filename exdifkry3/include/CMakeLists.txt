#From the build directory invoke using 
#cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
#cmake -DCMAKE_CXX_COMPILER=mpicxx ..
cmake_minimum_required(VERSION 3.10)
project(RadiativeTransfer2D LANGUAGES CXX)

# Use MPI compiler wrappers
# (Set externally via environment or command line)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Verbose output
set(CMAKE_VERBOSE_MAKEFILE ON)

# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# HYPRE
set(HYPRE_DIR "../../hypre/src/hypre")
include_directories(${HYPRE_DIR}/include)
link_directories(${HYPRE_DIR}/lib)

# Source files
file(GLOB SRC_FILES src/*.cpp)

# Executable
add_executable(radiative_transfer ${SRC_FILES})

# Link libraries
target_link_libraries(radiative_transfer ${MPI_CXX_LIBRARIES} HYPRE)

# Show paths
message(STATUS "MPI Include Path: ${MPI_INCLUDE_PATH}")
message(STATUS "HYPRE Include Path: ${HYPRE_DIR}/include")
message(STATUS "HYPRE Library Path: ${HYPRE_DIR}/lib")