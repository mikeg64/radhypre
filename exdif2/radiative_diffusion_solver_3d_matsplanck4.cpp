#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "mpi.h"

// Constants
const int N = 100;
const double L = 1.0;
const double dx = L / N;
const double dt = 1e-3;
const int max_time_steps = 100;
const int num_frequencies = 5;

// Function to write VTK output
void write_vtk(const std::string &filename, const std::vector<double> &temperature) {
    std::ofstream vtkFile(filename);
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Temperature field output\n";
    vtkFile << "ASCII\n";
    vtkFile << "DATASET STRUCTURED_POINTS\n";
    vtkFile << "DIMENSIONS " << N << " " << N << " " << N << "\n";
    vtkFile << "ORIGIN 0 0 0\n";
    vtkFile << "SPACING " << dx << " " << dx << " " << dx << "\n";
    vtkFile << "POINT_DATA " << N * N * N << "\n";
    vtkFile << "SCALARS Temperature float\n";
    vtkFile << "LOOKUP_TABLE default\n";

    for (double value : temperature) {
        vtkFile << value << "\n";
    }

    vtkFile.close();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int myid, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Initialize HYPRE structures
    HYPRE_IJMatrix A;
    HYPRE_IJVector rhs, solution;
    HYPRE_Solver solver;

    // Create matrix and vectors
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N*N*N-1, 0, N*N*N-1, &A);
    HYPRE_IJMatrixInitialize(A);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N*N*N-1, &rhs);
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N*N*N-1, &solution);
    HYPRE_IJVectorInitialize(rhs);
    HYPRE_IJVectorInitialize(solution);

    std::vector<double> temperature(N*N*N, 300.0); // Initial temperature

    for (int t = 0; t < max_time_steps; t++) {
        std::cout << "Time Step: " << t << std::endl;

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    int index = i*N*N + j*N + k;
                    double T_new = temperature[index] + (dt * 10.0 / 900.0); // Simplified heating term
                    temperature[index] = T_new;
                    HYPRE_IJVectorSetValues(rhs, 1, &index, &T_new);
                }
            }
        }

        // Solve system
        HYPRE_IJMatrixAssemble(A);
        HYPRE_IJVectorAssemble(rhs);
        HYPRE_IJVectorAssemble(solution);

        // Output VTK file at each step
        write_vtk("temperature_" + std::to_string(t) + ".vtk", temperature);
    }

    // Cleanup
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(rhs);
    HYPRE_IJVectorDestroy(solution);
    MPI_Finalize();
    
    return 0;
}