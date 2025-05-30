#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HYPRE_struct_ls.h"

const int N = 100;
const double L = 1.0;
const double dx = L / N;
const double h = 6.626e-34;  // Planck’s constant
const double c = 3.0e8;      // Speed of light
const double k_B = 1.38e-23; // Boltzmann constant

// Define Planckian source
double planck_source(double nu, double T) {
    return (2 * h * pow(nu, 3) / (c * c)) / (exp(h * nu / (k_B * T)) - 1);
}

// Define initial temperature distribution
double initial_temperature(int i, int j, int k) {
    return 300.0 + 50.0 * sin(i * dx) * cos(j * dx); // Example spatial variation
}

// Define specific heat capacity
double specific_heat_capacity(int i, int j, int k) {
    return 900.0 + 100.0 * exp(-k * dx); // Example material property variation
}

// Define absorption coefficient (spatially varying)
double absorption_coefficient(int i, int j, int k) {
    return 0.01 + 0.02 * cos(i * dx) * sin(j * dx);
}

int main() {
    // Initialize Hypre structures
    HYPRE_StructGrid grid;
    HYPRE_StructMatrix matrix;
    HYPRE_StructVector rhs, solution;
    HYPRE_StructSolver solver;

    // Grid setup
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);
    int ilower[3] = {0, 0, 0};
    int iupper[3] = {N - 1, N - 1, N - 1};
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);

    // Matrix setup
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, &matrix);
    HYPRE_StructMatrixInitialize(matrix);

    // Set matrix coefficients
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                double kappa = absorption_coefficient(i, j, k);
                double coeffs[7] = { -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, kappa };
                int stencil_indices[7] = { 0, 1, 2, 3, 4, 5, 6 };
                int index[3] = { i, j, k };
                HYPRE_StructMatrixSetValues(matrix, index, 7, stencil_indices, coeffs);
            }
        }
    }

    // Assemble matrix
    HYPRE_StructMatrixAssemble(matrix);

    // Create vectors
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &rhs);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &solution);
    HYPRE_StructVectorInitialize(rhs);
    HYPRE_StructVectorInitialize(solution);

    // Set RHS (Planckian source + heating update)
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                double T_old = initial_temperature(i, j, k);
                double cp = specific_heat_capacity(i, j, k);
                double source_intensity = planck_source(1.0e14, T_old);
                double energy_deposition = source_intensity * absorption_coefficient(i, j, k);
                double T_new = T_old + energy_deposition / cp;

                int index[3] = { i, j, k };
                HYPRE_StructVectorSetValues(rhs, index, &T_new);
            }
        }
    }

    // Assemble vectors
    HYPRE_StructVectorAssemble(rhs);
    HYPRE_StructVectorAssemble(solution);

    // Solve
    HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructPCGSetTol(solver, 1e-6);
    HYPRE_StructPCGSetup(solver, matrix, rhs, solution);
    HYPRE_StructPCGSolve(solver, matrix, rhs, solution);

    // Write results
    std::ofstream outFile("temperature_solution.txt");
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
    HYPRE_StructMatrixDestroy(matrix);
    HYPRE_StructGridDestroy(grid);

    return 0;
}