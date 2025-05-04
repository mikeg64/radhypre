
#include "ex.h"
#include <HYPRE.h>
#include <HYPRE_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_struct_ls.h>
#include <HYPRE_struct_mv.h>
//#include <mpi.h>
#include <math.h>
#include <stdio.h>

const double k_B = 1.380649e-23;  // Boltzmann constant
const double h = 6.62607015e-34;  // Planck constant
const double c = 3.0e8;           // Speed of light

double planck_law(double T, double lambda) {
    double numerator = 2.0 * h * pow(c, 2);
    double denominator = pow(lambda, 5) * (exp(h * c / (lambda * k_B * T)) - 1.0);
    return numerator / denominator;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Initialize HYPRE
    HYPRE_Init();

  
         printf("start \n");

    

    // Define problem size and coefficients
    int nx = 10, ny = 10, nz = 10;
    double heat_capacity = 1.0;        // Specific heat capacity
    double scattering_coeff = 0.01;    // Scattering coefficient
    double absorption_coeff = 0.1;     // Absorption coefficient
    double temperature = 300.0;        // Temperature in Kelvin
    double lambda = 500e-9;            // Wavelength in meters

    double ti = 0.0;                   // Initial time
    double tf = 1.0;                   // Final time
    double dt = 0.01;                  // Time step

    // Create grid
    HYPRE_StructGrid grid;

     printf("create grid \n");

    HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);

    // Define the grid dimensions
    HYPRE_Int ilower[3] = {0, 0, 0};
    HYPRE_Int iupper[3] = {nx-1, ny-1, nz-1};
    HYPRE_StructGridSetExtents(grid, ilower, iupper);

    // Define stencil
    HYPRE_StructStencil stencil;
    HYPRE_StructStencilCreate(3, 7, &stencil);
    // Define stencil entries (center, x, y, z directions)
    int offsets[7][3] = {{0, 0, 0}, {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    int offsetsf[21] = {0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1};
    for (int i = 0; i < 7; i++) {
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }

         printf("define stencil \n");

    // Create matrix and vector
    HYPRE_StructMatrix A;
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    HYPRE_StructVector b;
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
    HYPRE_StructVectorInitialize(b);

    HYPRE_StructVector x;
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    HYPRE_StructVectorInitialize(x);

    // Set up coefficients (example values)
    double h_grid = 1.0;
    double alpha = 1.0;
    double values[7] = {absorption_coeff + scattering_coeff + (1.0 / heat_capacity), 
                        -alpha / (h_grid * h_grid), -alpha / (h_grid * h_grid), 
                        -alpha / (h_grid * h_grid), -alpha / (h_grid * h_grid), 
                        -alpha / (h_grid * h_grid), -alpha / (h_grid * h_grid)};

    printf("time stepping loop \n");
    // Time-stepping loop
    for (double t = ti; t <= tf; t += dt) {
        // Loop over grid points to set matrix coefficients and RHS
        for (int i = ilower[0]; i <= iupper[0]; i++) {
            for (int j = ilower[1]; j <= iupper[1]; j++) {
                for (int k = ilower[2]; k <= iupper[2]; k++) {
                    HYPRE_Int idx[3] = {i, j, k};
                    HYPRE_StructMatrixSetValues(A, idx, 7, (HYPRE_Real*)offsets, values);
                    double b_val = (i == 0) ? planck_law(temperature, lambda) : 0.0; // Source term based on Planck law
                    HYPRE_StructVectorSetValues(b, idx, b_val);
                }
            }
        }

       printf("loop completed %f \n",t);
        // Finalize the matrix and vectors
        HYPRE_StructMatrixAssemble(A);
        HYPRE_StructVectorAssemble(b);
        HYPRE_StructVectorAssemble(x);

        printf("assembled\n");

        // Create solver
        HYPRE_StructSolver solver;
        HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructPCGSetMaxIter(solver, 100);
        HYPRE_StructPCGSetTol(solver, 1e-6);
        HYPRE_StructPCGSetTwoNorm(solver, 1);
        HYPRE_StructPCGSetPrintLevel(solver, 2);


        printf("created solver\n");

        // Setup and solve
        HYPRE_StructPCGSetup(solver, A, b, x);
        HYPRE_StructPCGSolve(solver, A, b, x);
        printf("solved\n");

        // Cleanup solver for the current time step
        HYPRE_StructPCGDestroy(solver);
    }

    // Cleanup
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructGridDestroy(grid);

    // Finalize HYPRE
    HYPRE_Finalize();
    MPI_Finalize();

    return 0;
}
