#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "HYPRE_struct_ls.h"
#include "HYPRE.h"
#include "mpi.h"

#define NX 50  // Grid points in X
#define NY 50  // Grid points in Y
#define TIME_STEPS 10 // Number of time steps
#define WAVELENGTHS 3 // Spectral bands
#define ALPHA 0.1 // Diffusion coefficient
#define h 6.626e-34 // Planck's constant
#define c 3.0e8 // Speed of light
#define k_B 1.38e-23 // Boltzmann constant

// Function to compute Planck’s spectral radiance
double planck_radiance(double T, double lambda) {
    return (2.0 * h * pow(c, 2)) / (pow(lambda, 5)) * (1.0 / (exp((h * c) / (lambda * k_B * T)) - 1.0));
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    HYPRE_Solver solver;

    int ilower[2] = {0, 0};
    int iupper[2] = {NX-1, NY-1};

    // Temperature field (initialized at room temperature)
    double T[NX][NY];
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            T[i][j] = 300.0;

    // Define wavelength bands (IR range)
    double lambda[WAVELENGTHS] = {1e-6, 5e-6, 10e-6};

    // Create grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);

    // Create stencil
    int stencil_size = 5;
    HYPRE_StructStencilCreate(2, stencil_size, &stencil);
    int offsets[5][2] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};

    for (int i = 0; i < stencil_size; i++)
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);

    // Create matrix
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, &A);
    HYPRE_StructMatrixSetStencil(A, stencil);
    HYPRE_StructMatrixInitialize(A);

    // Set matrix values
    double h_grid = 0.05;
    double values[5] = {1.0, -ALPHA / (h_grid * h_grid), -ALPHA / (h_grid * h_grid),
                        -ALPHA / (h_grid * h_grid), -ALPHA / (h_grid * h_grid)};

    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            int index[2] = {i, j};
            HYPRE_StructMatrixSetValues(A, index, stencil_size, values);
        }

    HYPRE_StructMatrixAssemble(A);

    // Create vectors
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);

    // Time loop
    for (int t = 0; t < TIME_STEPS; t++) {
        printf("Solving for time step %d...\n", t);

        // Update source using the **new temperature field**
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++) {
                int index[2] = {i, j};
                double source_val = 0.0;

                // Use the **temperature from the previous timestep**
                for (int w = 0; w < WAVELENGTHS; w++) {
                    source_val += planck_radiance(T[i][j], lambda[w]);
                }

                HYPRE_StructVectorSetValues(b, index, &source_val);
            }

        HYPRE_StructVectorAssemble(b);

        // Solve using PCG
        HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructPCGSetTol(solver, 1e-6);
        HYPRE_StructPCGSetMaxIter(solver, 100);
        HYPRE_StructPCGSetup(solver, A, b, x);
        HYPRE_StructPCGSolve(solver, A, b, x);

        // Update **temperature field** using the new solution with **spatial diffusion**
        double T_new[NX][NY];

        for (int i = 1; i < NX-1; i++)
            for (int j = 1; j < NY-1; j++) {
                int index[2] = {i, j};
                double sol;
                HYPRE_StructVectorGetValues(x, index, &sol);

                // Explicit diffusion update
                double diffusion = ALPHA * (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4 * T[i][j]);

                // Update temperature using both solution & diffusion
                T_new[i][j] = T[i][j] + sol * 0.01 + diffusion * 0.01;
            }

        // Copy new values into T
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                T[i][j] = T_new[i][j];

        // Destroy solver
        HYPRE_StructPCGDestroy(solver);
    }

    // Cleanup
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);

    MPI_Finalize();
    return 0;
}