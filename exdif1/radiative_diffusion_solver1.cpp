#include "ex.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include "HYPRE.h"
#include "HYPRE_struct_ls.h"

#define N 100   // Number of spatial points
#define L 1.0   // Length of the domain
#define PI 3.141592653589793

// Function to define a Gaussian source
double gaussian_source(int i, int N, double sigma) {
    double x = (double)i / (N - 1);
    return exp(-pow((x - 0.05) / sigma, 2)); // Centered near x=0.05
}

// Function to define spatially varying diffusion coefficient D(x)
double D_x(double x) {
    return 0.01 + 0.005 * sin(2 * PI * x); // Example variation
}

// Function to define spatially varying absorption coefficient kappa(x)
double kappa_x(double x) {
    return 0.02 + 0.01 * cos(2 * PI * x); // Example variation
}

int main() {
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    HYPRE_StructSolver solver;
    int ilower = 0, iupper = N - 1;

    // Initialize Hypre
    HYPRE_Init();

    // Create grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 1, &grid);
    HYPRE_StructGridSetExtents(grid, &ilower, &iupper);
    HYPRE_StructGridAssemble(grid);

    // Define stencil (5-point in 1D)
    HYPRE_StructStencilCreate(1, 3, &stencil);
    int offsets[3][1] = {{-1}, {0}, {1}};
    for (int i = 0; i < 3; i++) {
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }

    // Create matrix
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    // Set up coefficients with spatially varying D(x) and kappa(x)
    double dx = L / (N - 1);
    for (int i = 0; i < N; i++) {
        double x = (double)i / (N - 1);
        double D = D_x(x);
        double kappa = kappa_x(x);

        double coeffs[3] = {-D / (dx * dx), kappa + (2 * D / (dx * dx)), -D / (dx * dx)};
        int index[1] = {i};
        HYPRE_StructMatrixSetValues(A, index, 3, (int*)offsets, coeffs);
    }

    // Apply reflective boundary conditions at side walls
    int left_index[1] = {0};
    double left_coeffs[3] = {0, 1, -1};  // Zero flux at left boundary
    HYPRE_StructMatrixSetValues(A, left_index, 3, (int*)offsets, left_coeffs);

    // Apply open boundary at the end (Neumann condition)
    int right_index[1] = {N-1};
    double right_coeffs[3] = {-1, 1, 0}; // Allows radiation to leave
    HYPRE_StructMatrixSetValues(A, right_index, 3, (int*)offsets, right_coeffs);

    HYPRE_StructMatrixAssemble(A);

    // Create vectors
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);

    // Set Gaussian source with initial direction of 25 degrees
    double sigma = 0.02;
    for (int i = 0; i < N; i++) {
        double value = gaussian_source(i, N, sigma);
        int index[1] = {i};
        HYPRE_StructVectorSetValues(b, index, &value);
    }

    HYPRE_StructVectorAssemble(b);
    HYPRE_StructVectorAssemble(x);

    // Create solver (CG)
    HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_StructPCGSetMaxIter(solver, 1000);
    HYPRE_StructPCGSetTol(solver, 1e-6);
    HYPRE_StructPCGSetup(solver, A, b, x);
    HYPRE_StructPCGSolve(solver, A, b, x);

    // Write solution to file
    std::ofstream outfile("radiation_output.txt");
    double result[N];
    for (int i = 0; i < N; i++) {
        int index[1] = {i};
        HYPRE_StructVectorGetValues(x, index, &result[i]);
        outfile << i * dx << " " << result[i] << std::endl;
    }
    outfile.close();

    std::cout << "Solution saved to radiation_output.txt" << std::endl;

    // Clean up
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
    HYPRE_StructPCGDestroy(solver);
    HYPRE_Finalize();

    return 0;
}