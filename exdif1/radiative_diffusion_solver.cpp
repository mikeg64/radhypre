#include "ex.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <HYPRE_struct_ls.h>
#include <HYPRE_krylov.h>
#include <HYPRE.h>
#include <HYPRE_struct_mv.h>
#include <HYPRE_utilities.h>

// Define the problem parameters
const int N = 100; // Number of grid points
const double L = 1.0; // Length of the domain
const double dx = L / (N - 1); // Grid spacing
const double pi = 3.14159265358979323846;

// Function to compute the diffusion coefficient D(x)
double D(double x) {
    return 0.1 + 0.05 * sin(2 * pi * x);
}

// Function to compute the absorption coefficient kappa(x)
double kappa(double x) {
    return 1.0 + 0.5 * cos(2 * pi * x);
}

// Function to compute the source term S(x)
double S(double x) {
    double sigma = 0.01; // Width of the Gaussian
    double x0 = 0.1; // Center of the Gaussian
    return exp(-pow((x - x0) / sigma, 2));
}

int main(int argc, char *argv[]) {
    // Initialize HYPRE
    HYPRE_Init(argc, argv);

    // Create the grid
    HYPRE_StructGrid grid;
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 1, &grid);
    HYPRE_StructGridSetExtents(grid, {0}, {N-1});
    HYPRE_StructGridAssemble(grid);

    // Define the stencil
    HYPRE_StructStencil stencil;
    HYPRE_StructStencilCreate(1, 3, &stencil);
    HYPRE_StructStencilSetElement(stencil, 0, {0});
    HYPRE_StructStencilSetElement(stencil, 1, {-1});
    HYPRE_StructStencilSetElement(stencil, 2, {1});

    // Create the matrix
    HYPRE_StructMatrix A;
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    // Create the right-hand side vector
    HYPRE_StructVector b;
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
    HYPRE_StructVectorInitialize(b);

    // Create the solution vector
    HYPRE_StructVector x;
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    HYPRE_StructVectorInitialize(x);

    // Set up the matrix and vectors
    for (int i = 0; i < N; ++i) {
        double xi = i * dx;
        double Di = D(xi);
        double kappai = kappa(xi);
        double Si = S(xi);

        // Set the matrix values
        double values[3];
        int stencil_indices[3] = {0, 1, 2};
        values[0] = -2 * Di / (dx * dx) + kappai;
        values[1] = Di / (dx * dx);
        values[2] = Di / (dx * dx);
        HYPRE_StructMatrixSetValues(A, {i}, 3, stencil_indices, values);

        // Set the right-hand side values
        HYPRE_StructVectorSetValues(b, {i}, &Si);
    }

    // Apply reflective boundary conditions (Neumann)
    double zero = 0.0;
    HYPRE_StructMatrixSetValues(A, {0}, 1, {0}, &zero);
    HYPRE_StructMatrixSetValues(A, {N-1}, 1, {0}, &zero);
    HYPRE_StructVectorSetValues(b, {0}, &zero);
    HYPRE_StructVectorSetValues(b, {N-1}, &zero);

    // Assemble the matrix and vectors
    HYPRE_StructMatrixAssemble(A);
    HYPRE_StructVectorAssemble(b);
    HYPRE_StructVectorAssemble(x);

    // Set up and solve the linear system using BoomerAMG
    HYPRE_StructSolver solver;
    HYPRE_StructSolverCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructSolverSetTol(solver, 1e-6);
    HYPRE_StructSolverSetMaxIter(solver, 100);
    HYPRE_StructSolverSetPrintLevel(solver, 2);
    HYPRE_StructSolverSetup(solver, A, b, x);
    HYPRE_StructSolverSolve(solver, A, b, x);

    // Get the solution values
    std::vector<double> solution(N);
    for (int i = 0; i < N; ++i) {
        HYPRE_StructVectorGetValues(x, {i}, &solution[i]);
    }

    // Write the solution to a file
    std::ofstream outfile("solution.txt");
    for (int i = 0; i < N; ++i) {
        outfile << i * dx << " " << solution[i] << std::endl;
    }
    outfile.close();

    // Clean up
    HYPRE_StructSolverDestroy(solver);
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructGridDestroy(grid);
    HYPRE_Finalize();

    return 0;
}
