cmake_minimum_required(VERSION 3.14)
project(TypedMockExample LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Enable verbose output
set(CMAKE_VERBOSE_MAKEFILE ON)
# Find MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# HYPRE installation path
set(HYPRE_DIR "/home/mike/proj/hypre/src/hypre")
#set(GOOGLETEST_DIR "/home/mike/proj/googletest/build")
set(GOOGLETEST_DIR "/home/mike/proj/googletest")


#set(GOOGLETEST_DIR "/home/mike/proj/radhypre_classkrymgdif/exdifkry3/tests/googletest")
include_directories(${HYPRE_DIR}/include
                    ${GOOGLETEST_DIR}/googletest/include
                    ${GOOGLETEST_DIR}/googlemock/include)
link_directories(${HYPRE_DIR}/lib
                 ${GOOGLETEST_DIR}/build/lib)


add_subdirectory(${GOOGLETEST_DIR} ${CMAKE_BINARY_DIR}/gtest_build)


# Add GoogleTest + GoogleMock
#add_subdirectory(/home/mike/proj/googletest ${CMAKE_BINARY_DIR}/gtest_build)

# Test executable
add_executable(typed_mock_test
    typedmocktest.cpp
    vecwrapperadditiontest.cpp
)


target_include_directories(typed_mock_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(typed_mock_test
    gmock
    gmock_main
    gtest
    gtest_main
)

enable_testing()
add_test(NAME typed_mock_test COMMAND typed_mock_test)