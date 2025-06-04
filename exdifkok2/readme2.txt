Cmake tutorial
https://cliutils.gitlab.io/modern-cmake/chapters/basics/example.html


https://github.com/kokkos/kokkos-tutorials

https://github.com/kokkos/kokkos

https://kokkos.org/kokkos-core-wiki/


Building kokkos from kokkos distribution


The raw Makefiles require Makefile variables to be properly configured. In most examples, this is KOKKOS_PATH pointing to the Kokkos source directory and KOKKOS_DEVICES

export KOKKOS_PATH=pathtothe distributionbuiltfrommaster

make -j KOKKOS_DEVICES=OpenMP
make -j KOKKOS_DEVICES=Cuda

Olympia project folders built to cuda and host folder
Module load gcc/9.2.0
Module load cuda/11.4

To build using cmake the CMakeLists.txt must include



find_package(Kokkos REQUIRED)

target_link_libraries(hello_world PRIVATE Kokkos::kokkos)


Run cmake using

cmake . -DKokkos_ROOT=$KOKKOS_PATH


