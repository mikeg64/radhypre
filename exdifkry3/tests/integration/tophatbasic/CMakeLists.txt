
# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# HYPRE installation path
set(HYPRE_DIR "/home/mike/proj/hypre/src/hypre")
include_directories(${HYPRE_DIR}/include)
link_directories(${HYPRE_DIR}/lib)



add_executable(radt-crookedpipe main.cpp)
target_link_libraries(radt-crookedpipe  ${MPI_CXX_LIBRARIES} radt_lib HYPRE)
target_include_directories(radt-crookedpipe PUBLIC ${CMAKE_SOURCE_DIR}/include)