
#include "ex.h"
#include <HYPRE.h>
#include <HYPRE_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE_struct_ls.h>
#include <HYPRE_struct_mv.h>
//#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const double k_B = 1.380649e-23;  // Boltzmann constant
const double h = 6.62607015e-34;  // Planck constant
const double c = 3.0e8;           // Speed of light

double planck_law(double T, double lambda) {
    double numerator = 2.0 * h * pow(c, 2);
    double denominator = pow(lambda, 5) * (exp(h * c / (lambda * k_B * T)) - 1.0);
    return numerator / denominator;
}



int main(int argc, char *argv[])
{
    
    MPI_Init(&argc, &argv);

    // Initialize HYPRE
    HYPRE_Init();
    // Define matrix and grid dimensions
    HYPRE_StructMatrix A;
    HYPRE_StructGrid grid;
    
    // Create structured grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);

    // Set grid extents (example 3D case)
    int ilower[3] = {0, 0, 0};
    int iupper[3] = {9, 9, 9};  // Assuming a 10x10x10 grid
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    
    // Finalize grid setup
    HYPRE_StructGridAssemble(grid);


    // Set stencil size (for 7-point stencil)
    int stencil_size = 7;
    //HYPRE_StructStencil stencil;
    //HYPRE_StructStencilCreate(3, stencil_size, &stencil);

    // Define offsets for stencil
    int offsets[7][3] = {
        { 0,  0,  0},  // Center
        { 1,  0,  0},  // Right
        {-1,  0,  0},  // Left
        { 0,  1,  0},  // Front
        { 0, -1,  0},  // Back
        { 0,  0,  1},  // Up
        { 0,  0, -1}   // Down
    };

    // Define stencil
    HYPRE_StructStencil stencil;
    HYPRE_StructStencilCreate(3, 7, &stencil);
    // Define stencil entries (center, x, y, z directions)
    for (int i = 0; i < 7; i++) {
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }





    // Create structured matrix
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid,stencil, &A);

    // Assemble stencil
    //HYPRE_StructStencilAssemble(stencil);

    // Define coefficients
    double absorption_coeff = 0.5;
    double scattering_coeff = 0.2;
    double heat_capacity = 1.0;
    double alpha = 0.1;
    double h_grid = 0.05;

    // Define values for stencil
    double values[7] = {
        absorption_coeff + scattering_coeff + (1.0 / heat_capacity),
        -alpha / (h_grid * h_grid),
        -alpha / (h_grid * h_grid),
        -alpha / (h_grid * h_grid),
        -alpha / (h_grid * h_grid),
        -alpha / (h_grid * h_grid),
        -alpha / (h_grid * h_grid)
    };

    // Set matrix values in a loop
    for (int i = ilower[0]; i <= iupper[0]; i++) {
        for (int j = ilower[1]; j <= iupper[1]; j++) {
            for (int k = ilower[2]; k <= iupper[2]; k++) {
                int index[3] = {i, j, k};  // Current grid point
                //HYPRE_StructMatrixSetValues(A, index, stencil_size, values);
                HYPRE_StructMatrixSetValues(A, index, 7, (HYPRE_Real*)offsets, values);
            }
        }
    }

    // Assemble matrix
    HYPRE_StructMatrixAssemble(A);

    // Finalize
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);

    return 0;
}