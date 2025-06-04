A hybrid OpenMPI + Kokkos + HYPRE example combines distributed memory parallelism (via MPI), shared memory or accelerator parallelism (via Kokkos), and scalable linear solvers (via HYPRE). While there isn't a single public example that integrates all three in a minimal form, here's a simplified template to get you started:

🧩 Overview of the Example

MPI: Handles domain decomposition across nodes.
Kokkos: Manages on-node parallelism (e.g., GPU or multicore CPU).
HYPRE: Solves a linear system (e.g., from a discretized PDE).
📁 Project Structure

hybrid_solver/

├── CMakeLists.txt

└── src/

    └── main.cpp

📄 CMakeLists.txt

cmake_minimum_required(VERSION 3.10)

project(HybridSolver)

 

find_package(MPI REQUIRED)

find_package(Kokkos REQUIRED)

find_package(HYPRE REQUIRED)

 

add_executable(hybrid_solver src/main.cpp)

target_link_libraries(hybrid_solver MPI::MPI_CXX Kokkos::kokkos HYPRE::HYPRE)

📄 src/main.cpp

#include <mpi.h>

#include <Kokkos_Core.hpp>

#include "HYPRE.h"

#include "HYPRE_parcsr_ls.h"

 

int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    Kokkos::initialize(argc, argv);

 

    {

        int rank, size;

        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        MPI_Comm_size(MPI_COMM_WORLD, &size);

 

        if (rank == 0) {

            std::cout << "Running hybrid MPI + Kokkos + HYPRE example on " << size << " processes.\n";

        }

 

        // Kokkos parallel region

        Kokkos::parallel_for("HelloKokkos", 10, KOKKOS_LAMBDA(int i) {

            printf("Hello from Kokkos on rank %d, i = %d\\n", rank, i);

        });

 

        // HYPRE setup (simplified)

        HYPRE_IJMatrix A;

        HYPRE_IJVector b, x;

        HYPRE_ParCSRMatrix par_A;

        HYPRE_ParVector par_b, par_x;

 

        HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, 9, 0, 9, &A);

        HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

        HYPRE_IJMatrixInitialize(A);

        // Fill matrix A here...

 

        HYPRE_IJMatrixAssemble(A);

        HYPRE_IJMatrixGetObject(A, (void**) &par_A);

 

        // Solver setup

        HYPRE_Solver solver;

        HYPRE_BoomerAMGCreate(&solver);

        HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x);

        HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x);

 

        HYPRE_BoomerAMGDestroy(solver);

        HYPRE_IJMatrixDestroy(A);

    }

 

    Kokkos::finalize();

    MPI_Finalize();

    return 0;

}

Would you like me to package this into a downloadable zip file or extend it to solve a real PDE (e.g., Poisson equation) using domain decomposition?

 

 

 

 

 

 

To compile and run the hybrid Kokkos + MPI + HYPRE program, follow these steps. This assumes you have installed:

MPI (e.g., OpenMPI or MPICH)
Kokkos (with the desired backend like CUDA, HIP, or OpenMP)
HYPRE (built with MPI and optionally GPU support)
🛠️ Step 1: Create the Build Directory

mkdir build

cd build

🧱 Step 2: Configure with CMake

cmake .. \

  -DCMAKE_CXX_COMPILER=mpicxx \

  -DKokkos_DIR=/path/to/kokkos/install/lib/cmake/Kokkos \

  -DHYPRE_DIR=/path/to/hypre/install

Replace the paths with the actual installation directories of Kokkos and HYPRE.

🧪 Step 3: Build the Executable

make -j

This should produce an executable named hybrid_solver.

🚀 Step 4: Run with MPI

mpirun -np 4 ./hybrid_solver

This runs the program on 4 MPI processes. You can adjust the number of processes depending on your system.

⚙️ Optional: Enable GPU Backend

If you built Kokkos with CUDA or HIP, you can enable GPU execution by setting environment variables:

export OMP_NUM_THREADS=1  # if using OpenMP

export CUDA_VISIBLE_DEVICES=0  # if using CUDA

Would you like a script to automate this build and run process, or help with installing any of the dependencies?

 