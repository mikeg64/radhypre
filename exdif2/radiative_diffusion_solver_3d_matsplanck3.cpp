#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_parcsr_mv.h"
#include "mpi.h"

// Constants
const int N = 100;           // Grid resolution
const double L = 1.0;        // Box length
const double dx = L / N;     // Element size
const double dt = 1e-3;      // Time step size
const int max_time_steps = 100;
const int num_frequencies = 5;
const double h = 6.626e-34;  // Planck’s constant
const double c = 3.0e8;      // Speed of light
const double k_B = 1.38e-23; // Boltzmann constant

// Sample frequency list
double frequencies[num_frequencies] = { 1.0e14, 2.0e14, 3.0e14, 4.0e14, 5.0e14 };

// Planck source function
double planck_source(double nu, double T) {
    return (2 * h * pow(nu, 3) / (c * c)) / (exp(h * nu / (k_B * T)) - 1);
}

// Initial temperature distribution
double initial_temperature(int i, int j, int k) {
    return 300.0 + 50.0 * sin(i * dx) * cos(j * dx);
}

// Specific heat capacity variation
double specific_heat_capacity(int i, int j, int k) {
    return 900.0 + 100.0 * exp(-k * dx);
}

// Absorption coefficient variation
double absorption_coefficient(int i, int j, int k, double nu) {
    return 0.01 + 0.02 * cos(i * dx) * sin(j * dx) * exp(-nu / 1.0e14);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int myid, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Initialize HYPRE matrix solver
    HYPRE_IJMatrix A;
    HYPRE_IJVector rhs, solution;
    HYPRE_Solver solver;

    // Create matrix and vectors
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N*N*N-1, 0, N*N*N-1, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N*N*N-1, &rhs);
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N*N*N-1, &solution);
    HYPRE_IJVectorSetObjectType(rhs, HYPRE_PARCSR);
    HYPRE_IJVectorSetObjectType(solution, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(rhs);
    HYPRE_IJVectorInitialize(solution);

    // Time stepping loop
    for (int t = 0; t < max_time_steps; t++) {
        std::cout << "Time Step: " << t << std::endl;

        // Frequency loop
        for (int nu_idx = 0; nu_idx < num_frequencies; nu_idx++) {
            double nu = frequencies[nu_idx];

            // Assemble system
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < N; k++) {
                        int index = i*N*N + j*N + k;
                        double T_old = initial_temperature(i, j, k);
                        double cp = specific_heat_capacity(i, j, k);
                        double source_intensity = planck_source(nu, T_old);
                        double energy_deposition = source_intensity * absorption_coefficient(i, j, k, nu);
                        double T_new = T_old + (dt * energy_deposition / cp);

                        double values[7] = { -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 1.0 };
                        int stencil_indices[7] = { index-N*N, index-N, index-1, index+1, index+N, index+N*N, index };

                        HYPRE_IJMatrixSetValues(A, 1, &index, 7, stencil_indices, values);
                        HYPRE_IJVectorSetValues(rhs, 1, &index, &T_new);
                    }
                }
            }

            // Solve system
            HYPRE_IJMatrixAssemble(A);
            HYPRE_IJVectorAssemble(rhs);
            HYPRE_IJVectorAssemble(solution);

            HYPRE_ParCSRMatrix par_A;
            HYPRE_ParCSRVector par_rhs, par_sol;
            HYPRE_IJMatrixGetObject(A, (void **)&par_A);
            HYPRE_IJVectorGetObject(rhs, (void **)&par_rhs);
            HYPRE_IJVectorGetObject(solution, (void **)&par_sol);

            HYPRE_BoomerAMGCreate(&solver);
            HYPRE_BoomerAMGSetTol(solver, 1e-6);
            HYPRE_BoomerAMGSetup(solver, par_A, par_rhs, par_sol);
            HYPRE_BoomerAMGSolve(solver, par_A, par_rhs, par_sol);

            // Write solution
            std::ofstream outFile("hypre_fem_solution_" + std::to_string(t) + "_" + std::to_string(nu_idx) + ".txt");
            for (int i = 0; i < N*N*N; i++) {
                double value;
                HYPRE_IJVectorGetValues(solution, 1, &i, &value);
                outFile << i << " " << value << std::endl;
            }
            outFile.close();

            // Reset matrix for next step
            HYPRE_IJMatrixDestroy(A);
            HYPRE_IJVectorDestroy(rhs);
            HYPRE_IJVectorDestroy(solution);
        }
    }

    // Cleanup
    HYPRE_BoomerAMGDestroy(solver);
    MPI_Finalize();
    return 0;
}