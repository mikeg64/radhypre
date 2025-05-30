#include "ex.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>
#include <HYPRE_struct_ls.h>
#include <HYPRE_struct_mv.h>
#include <mpi.h>

// Constants
const int N = 100; // Grid size
const double L = 1.0; // Domain length
const double dx = L / (N - 1); // Grid spacing
const double dt = 1e-4; // Time step
const double t_end = 1.0; // End time
const int num_freq_bins = 10; // Number of frequency bins
const double c = 3e10; // Speed of light in cm/s
const double a = 7.5646e-15; // Radiation constant in erg/cm^3/K^4
const double sigma = 5.6704e-5; // Stefan-Boltzmann constant in erg/cm^2/s/K^4

// Function to compute Planck distribution
double B_nu(double T, double nu) {
    return (2 * h * pow(nu, 3) / pow(c, 2)) / (exp(h * nu / (k * T)) - 1);
}

// Function to compute diffusion coefficient
double D(double x, double y, double z) {
    return 0.1 + 0.05 * sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
}

// Function to compute absorption coefficient
double kappa(double x, double y, double z) {
    return 1.0 + 0.5 * cos(2 * M_PI * x) * cos(2 * M_PI * y) * cos(2 * M_PI * z);
}

// Function to initialize temperature field
void initialize_temperature(std::vector<std::vector<std::vector<double>>> &T) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                T[i][j][k] = 300.0; // Initial temperature in K
            }
        }
    }
}

// Function to initialize radiation field
void initialize_radiation(std::vector<std::vector<std::vector<std::vector<double>>>> &phi, const std::vector<std::vector<std::vector<double>>> &T) {
    for (int n = 0; n < num_freq_bins; ++n) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    phi[n][i][j][k] = a * pow(T[i][j][k], 4); // Initial radiation energy density
                }
            }
        }
    }
}

// Function to solve the radiative diffusion equation for a single frequency bin using HYPRE
void solveRadiativeDiffusion(int n, std::vector<std::vector<std::vector<double>>> &phi_n, const std::vector<std::vector<std::vector<double>>> &T) {
    // HYPRE setup
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    HYPRE_StructSolver solver;
    HYPRE_StructSolver precond;

    // Create an empty 3D grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);

    // Define the grid size
    int ilower[3] = {0, 0, 0};
    int iupper[3] = {N-1, N-1, N-1};
    HYPRE_StructGridSetExtents(grid, ilower, iupper);

    // Add periodic boundary conditions
    int periodic[3] = {0, 0, 0};
    HYPRE_StructGridSetPeriodic(grid, periodic);

    // Assemble the grid
    HYPRE_StructGridAssemble(grid);

    // Create a 7-point stencil
    HYPRE_StructStencilCreate(3, 7, &stencil);
    int offsets[7][3] = {{0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    for (int i = 0; i < 7; ++i) {
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }

    // Create the matrix
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    // Create the right-hand side and solution vectors
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);

    // Set matrix coefficients
    double values[7];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                int index[3] = {i, j, k};
                double D_val = D(i * dx, j * dx, k * dx);
                double kappa_val = kappa(i * dx, j * dx, k * dx);
                values[0] = -6.0 * D_val / (dx * dx) + kappa_val;
                values[1] = values[2] = values[3] = values[4] = values[5] = values[6] = D_val / (dx * dx);
                HYPRE_StructMatrixSetValues(A, index, 7, offsets, values);
            }
        }
    }

    // Assemble the matrix
    HYPRE_StructMatrixAssemble(A);

    // Set right-hand side values
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                int index[3] = {i, j, k};
                double kappa_val = kappa(i * dx, j * dx, k * dx);
                double B_val = B_nu(T[i][j][k], n);
                double rhs_val = kappa_val * B_val;
                HYPRE_StructVectorSetValues(b, index, &rhs_val);
            }
        }
    }

    // Assemble the right-hand side vector
    HYPRE_StructVectorAssemble(b);

    // Set initial guess for the solution
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                int index[3] = {i, j, k};
                double phi_val = phi_n[i][j][k];
                HYPRE_StructVectorSetValues(x, index, &phi_val);
            }
        }
    }

    // Assemble the solution vector
    HYPRE_StructVectorAssemble(x);

    // Create and set up the solver
    HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructPCGSetTol(solver, 1e-6);
    HYPRE_StructPCGSetMaxIter(solver, 100);
    HYPRE_StructPCGSetTwoNorm(solver, 1);
    HYPRE_StructPCGSetPrintLevel(solver, 2);
    HYPRE_StructPCGSetLogging(solver, 1);

    // Create and set up the preconditioner
    HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
    HYPRE_StructSMGSetMemoryUse(precond, 0);
    HYPRE_StructSMGSetMaxIter(precond, 1);
    HYPRE_StructSMGSetTol(precond, 0.0);
    HYPRE_StructSMGSetZeroGuess(precond);
    HYPRE_StructSMGSetNumPreRelax(precond, 1);
    HYPRE_StructSMGSetNumPostRelax(precond, 1);

    // Set the preconditioner
    HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve, HYPRE_StructSMGSetup, precond);

    // Set up and solve the system
    HYPRE_StructPCGSetup(solver, A, b, x);
    HYPRE_StructPCGSolve(solver, A, b, x);

    // Get the solution values
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                int index[3] = {i, j, k};
                double phi_val;
                HYPRE_StructVectorGetValues(x, index, &phi_val);
                phi_n[i][j][k] = phi_val;
            }
        }
    }

    // Clean up
    HYPRE_StructPCGDestroy(solver);
    HYPRE_StructSMGDestroy(precond);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructGridDestroy(grid);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    HYPRE_Init();

    // Initialize temperature and radiation fields
    std::vector<std::vector<std::vector<double>>> T(N, std::vector<std::vector<double>>(N, std::vector<double>(N)));
    std::vector<std::vector<std::vector<std::vector<double>>>> phi(num_freq_bins, std::vector<std::vector<std::vector<double>>>(N, std::vector<std::vector<double>>(N, std::vector<double>(N))));

    initialize_temperature(T);
    initialize_radiation(phi, T);

    // Time loop
    for (double t = 0; t < t_end; t += dt) {
        // Nonlinear iteration for coupling
        for (int iter = 0; iter < 10; ++iter) {
            // Solve radiative diffusion equation for each frequency bin
            for (int n = 0; n < num_freq_bins; ++n) {
                solveRadiativeDiffusion(n, phi[n], T);
            }

            // Update temperature field
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    for (int k = 0; k < N; ++k) {
                        double phi_sum = 0.0;
                        for (int n = 0; n < num_freq_bins; ++n) {
                            phi_sum += phi[n][i][j][k];
                        }
                        T[i][j][k] += dt * (phi_sum - a * pow(T[i][j][k], 4)) / (rho * c_v);
                    }
                }
            }
        }

        // Output temperature and radiation fields
        std::ofstream temp_file("temperature_" + std::to_string(t) + ".txt");
        std::ofstream rad_file("radiation_" + std::to_string(t) + ".txt");

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    temp_file << T[i][j][k] << " ";
                    for (int n = 0; n < num_freq_bins; ++n) {
                        rad_file << phi[n][i][j][k] << " ";
                    }
                    rad_file << std::endl;
                }
                temp_file << std::endl;
            }
            temp_file << std::endl;
        }

        temp_file.close();
        rad_file.close();
    }

    HYPRE_Finalize();
    MPI_Finalize();

    return 0;
}
