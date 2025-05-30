#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HYPRE_struct_ls.h"

// Constants
const int N = 100;
const double L = 1.0;
const double dx = L / N;
const double dt = 1e-3; // Time step size
const int max_time_steps = 100; // Number of time steps

const int num_frequencies = 5; // Simulate over 5 radiation frequencies
const double h = 6.626e-34;  // Planck’s constant
const double c = 3.0e8;      // Speed of light
const double k_B = 1.38e-23; // Boltzmann constant

// Frequency list (sample)
double frequencies[num_frequencies] = { 1.0e14, 2.0e14, 3.0e14, 4.0e14, 5.0e14 };

// Planckian source for each frequency
double planck_source(double nu, double T) {
    return (2 * h * pow(nu, 3) / (c * c)) / (exp(h * nu / (k_B * T)) - 1);
}

// Initial temperature distribution
double initial_temperature(int i, int j, int k) {
    return 300.0 + 50.0 * sin(i * dx) * cos(j * dx); // Example variation
}

// Specific heat capacity variation
double specific_heat_capacity(int i, int j, int k) {
    return 900.0 + 100.0 * exp(-k * dx); // Example material property variation
}

// Absorption coefficient (spatially varying)
double absorption_coefficient(int i, int j, int k, double nu) {
    return 0.01 + 0.02 * cos(i * dx) * sin(j * dx) * exp(-nu / 1.0e14);
}

int main() {
    // Initialize Hypre
    HYPRE_StructGrid grid;
    HYPRE_StructMatrix matrix;
    HYPRE_StructVector rhs, solution;
    HYPRE_StructSolver solver;

    // Create grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);
    int ilower[3] = {0, 0, 0};
    int iupper[3] = {N - 1, N - 1, N - 1};
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);

    // Create matrix
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, &matrix);
    HYPRE_StructMatrixInitialize(matrix);

    // Simulation over multiple time steps
    for (int t = 0; t < max_time_steps; t++) {
        std::cout << "Time Step: " << t << std::endl;

        // Update temperature and radiation at each time step
        for (int nu_idx = 0; nu_idx < num_frequencies; nu_idx++) {
            double nu = frequencies[nu_idx];

            // Create vectors
            HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhs);
            HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solution);
            HYPRE_StructVectorInitialize(rhs);
            HYPRE_StructVectorInitialize(solution);

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < N; k++) {
                        double T_old = initial_temperature(i, j, k);
                        double cp = specific_heat_capacity(i, j, k);
                        double source_intensity = planck_source(nu, T_old);
                        double energy_deposition = source_intensity * absorption_coefficient(i, j, k, nu);
                        double T_new = T_old + (dt * energy_deposition / cp);

                        int index[3] = { i, j, k };
                        HYPRE_StructVectorSetValues(rhs, index, &T_new);
                    }
                }
            }

            // Solve the system
            HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
            HYPRE_StructPCGSetTol(solver, 1e-6);
            HYPRE_StructPCGSetup(solver, matrix, rhs, solution);
            HYPRE_StructPCGSolve(solver, matrix, rhs, solution);

            // Save results for this frequency and time step
            std::ofstream outFile("solution_" + std::to_string(t) + "_" + std::to_string(nu_idx) + ".txt");
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) {
                    for (int k = 0; k < N; k++) {
                        int index[3] = { i, j, k };
                        double value;
                        HYPRE_StructVectorGetValues(solution, index, &value);
                        outFile << i << " " << j << " " << k << " " << value << std::endl;
                    }
                }
            }
            outFile.close();

            // Cleanup
            HYPRE_StructPCGDestroy(solver);
            HYPRE_StructVectorDestroy(rhs);
            HYPRE_StructVectorDestroy(solution);
        }
    }

    // Final Cleanup
    HYPRE_StructMatrixDestroy(matrix);
    HYPRE_StructGridDestroy(grid);

    return 0;
}