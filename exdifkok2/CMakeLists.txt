cmake_minimum_required(VERSION 3.10)

project(HybridSolver)

 

find_package(MPI REQUIRED)

find_package(Kokkos REQUIRED)

find_package(HYPRE REQUIRED)

 

add_executable(hybrid_solver src/main.cpp)

target_link_libraries(hybrid_solver MPI::MPI_CXX Kokkos::kokkos HYPRE::HYPRE)