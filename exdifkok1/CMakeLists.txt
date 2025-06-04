cmake_minimum_required(VERSION 3.10)

project(DiffusiveSolver)

 

# Find Kokkos

find_package(Kokkos REQUIRED)

 

# Find HYPRE

find_package(HYPRE REQUIRED)

 

add_executable(diffusive_solver src/main.cpp)

 

# Link Kokkos and HYPRE

target_link_libraries(diffusive_solver Kokkos::kokkos HYPRE::HYPRE)