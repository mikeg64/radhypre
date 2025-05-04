

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ex.h"
#include "HYPRE_struct_ls.h"
#include "HYPRE.h"
#include "mpi.h"

#define NX 50  // Grid points in X
#define NY 50  // Grid points in Y
#define TIME_STEPS 100000 // Number of time steps
#define WAVELENGTHS 3 // Spectral bands
#define ALPHA 0.2 // Diffusion coefficient
#define h 6.626e-34 // Planck's constant
#define c 3.0e8 // Speed of light
#define k_B 1.38e-23 // Boltzmann constant
#define DELTAX 0.001

// Radiation entrance settings
#define APERTURE_START NY/3  // Aperture position start
#define APERTURE_END 2*NY/3  // Aperture position end

// Fixed angle for directional radiation (Adjustable parameter)
#define ANGLE_DEGREES 90.0  // Change this to vary direction

// Convert angle to bias factor (favor propagation direction)
double compute_biasx_factor(double angle) {
    return fabs(cos(angle * M_PI / 180.0));  // Normalize bias effect
    printf("m_pik %f\n",M_PI);
}

double compute_biasy_factor(double angle) {
    return fabs(sin(angle * M_PI / 180.0));  // Normalize bias effect
}

// Function to compute Planck’s spectral radiance
double planck_radiance(double T, double lambda) {
    return (2.0 * h * pow(c, 2)) / (pow(lambda, 5)) * (1.0 / (exp((h * c) / (lambda * k_B * T)) - 1.0));
}




void write_vtk_file(double T[NX][NY], int timestep) {
    FILE *vtkFile;
    char filename[50];

    sprintf(filename, "radiation_output_%d.vtk", timestep);
    vtkFile = fopen(filename, "w");

    if (!vtkFile) {
        printf("Error: Could not open file for writing!\n");
        return;
    }

    // Write VTK header
    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Radiation Transport Simulation\n");
    fprintf(vtkFile, "ASCII\n");
    fprintf(vtkFile, "DATASET STRUCTURED_POINTS\n");
    fprintf(vtkFile, "DIMENSIONS %d %d 1\n", NX, NY);
    fprintf(vtkFile, "ORIGIN 0 0 0\n");
    fprintf(vtkFile, "SPACING 1 1 1\n");

    // Write data section
    fprintf(vtkFile, "POINT_DATA %d\n", (NX * NY));
    fprintf(vtkFile, "SCALARS Temperature float 1\n");
    fprintf(vtkFile, "LOOKUP_TABLE default\n");
 
    // Write temperature values in row-major order
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            fprintf(vtkFile, "%f ", T[i][j]);
        }
        fprintf(vtkFile, "\n");
    }

    fclose(vtkFile);
    printf("VTK file written successfully: %s\n", filename);
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


    for (int j = (NY/2)-3; j < (NY/2)+3; j++)
    // for (int j = 3; j < 6; j++)
            for (int i = 0; i < 3; i++)
                     T[i][j] = 50000.0; // Initial temperature (room temp)    

    // Define wavelength bands (IR range)
    double lambda[WAVELENGTHS] = {1e-6, 5e-6, 10e-6};

    // Compute bias factor based on angle
    double bias_x = 1.2*compute_biasx_factor(ANGLE_DEGREES);
    double bias_y = 1.2*compute_biasy_factor(ANGLE_DEGREES);

    // Create grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);

    // Create stencil
    int stencil_size = 5;
    HYPRE_StructStencilCreate(2, stencil_size, &stencil);
    int offsets[5][2] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};  // 5-point diffusion stencil

    for (int i = 0; i < stencil_size; i++)
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);

    // Create matrix
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
    //HYPRE_StructMatrixSetStencil(A, stencil);
    HYPRE_StructMatrixInitialize(A);

    // Set matrix values for diffusion-based transport
    double h_grid = 0.05;
    double values[5] = {1.0, -ALPHA*bias_x / (h_grid * h_grid), -ALPHA / (h_grid * h_grid),
                        -ALPHA*bias_y / (h_grid * h_grid), -ALPHA / (h_grid * h_grid)};

    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            int index[2] = {i, j};
            //HYPRE_StructMatrixSetValues(A, index, stencil_size, values);
            HYPRE_StructMatrixSetValues(A, index, 7, (HYPRE_Real*)offsets, values);
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

                // Only emit radiation if within the aperture region at x = 0
                if (i == 0 && j >= APERTURE_START && j <= APERTURE_END) {
                    for (int w = 0; w < WAVELENGTHS; w++) {
                        source_val += planck_radiance(T[i][j], lambda[w]);
                    }
                }

                HYPRE_StructVectorSetValues(b, index, source_val);
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

        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++) {
                int index[2] = {i, j};
                double sol;
                HYPRE_StructVectorGetValues(x, index, &sol);

                // Directional bias update based on user-defined angle
                //double directional_term = bias_factor * (T[i+1][j+1] - T[i][j]); // Favor movement at angle

                //double diffusion = ALPHA * (T[i+1][j] + T[i-1][j] + T[i][j+1] + T[i][j-1] - 4 * T[i][j]);

                // Final temperature update
                T_new[i][j] = T[i][j] + sol * 0.025 ;
            }


        for (int i = 1; i < NX-1; i++)
            for (int j = 1; j < NY-1; j++) {
                int index[2] = {i, j};
                double sol;
                HYPRE_StructVectorGetValues(x, index, &sol);

                // Directional bias update based on user-defined angle
                double directional_term = 0.0; // -bias_factor * (T[i+1][j+1] - T[i][j]); // Favor movement at angle

                double diffusion = ALPHA * (T[i+1][j] + 1.1*bias_x*T[i-1][j] + T[i][j+1] +1.1*bias_y* T[i][j-1] - 4 * T[i][j]);

                // Final temperature update
                T_new[i][j] = T[i][j] + diffusion * 0.025 + directional_term;
            }


        // Implement **reflective boundary conditions**
        for (int j = 0; j < NY; j++) {
            T_new[0][j] = T[0][j]; // Reflective at entrance
            T_new[NX-1][j] *= 0.9; // Exit wall open (allows energy escape)
        }


        // Implement **reflective boundary conditions** for side walls
        for (int i = 0; i < NX; i++) {
            T_new[i][0] = T[i][1];      // Reflect left boundary (j=0)
            T_new[i][NY-1] = T[i][NY-2]; // Reflect right boundary (j=NY-1)
        }

        //for (int j = 3; j < 6; j++)
        for (int j = (NY/2)-3; j < (NY/2)+3; j++)
            for (int i = 0; i < 3; i++)
                     T_new[i][j] = 50000.0; // Initial temperature (room temp)    

        
        // Copy new values into T
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
            {
                T[i][j] = T_new[i][j];
                if(t%100==0)
                    printf("Timestep %d: T[%d, %d] = %lf\n", t, i, j, T[i][j]);
            }

                // Write VTK output for visualization
        if(t%100==0)
            write_vtk_file(T, t/100);


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