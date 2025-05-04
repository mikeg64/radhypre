#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "HYPRE.h"
#include "HYPRE_struct_ls.h"

#define NX 50   // Grid points in X
#define NY 50   // Grid points in Y
#define ANGLES 8  // Discrete ordinates directions
#define TIME_STEPS 1000  // Number of time steps
#define DT 0.0000000001  // Time step size
#define ALPHA 0.00001  // Absorption coefficient
#define h 6.626e-34  
#define c 3.0e8  
#define k_B 1.38e-23  

// Radiation entrance settings
#define APERTURE_START NY/6  // Aperture position start
#define APERTURE_END 2*NY/6  // Aperture position end



// Define discrete angles (SN method)
double theta[ANGLES] = {0, 45, 90, 135, 180, 225, 270, 315};

// Radiation intensity field per angle
double I[NX][NY][ANGLES];

double heat_capacity(double T)  {
    return 0.0001;  // Placeholder for specific heat capacity
}

double sigmaa(double T)  {
    return 0.01;  // Placeholder for absorption
}

// Function to compute Planck’s spectral radiance
double planck_radiance(double T) {
    return (2.0 * h * pow(c, 2)) / (pow(1e-6, 5)) * (1.0 / (exp((h * c) / (1e-6 * k_B * T)) - 1.0));
}

// Function to compute Planck’s spectral radiance
//double planck_radiance(double T, double lambda) {
//    return (2.0 * h * pow(c, 2)) / (pow(lambda, 5)) * (1.0 / (exp((h * c) / (lambda * k_B * T)) - 1.0));
//}


void introduce_source_term(double I[NX][NY][ANGLES]) {
    double T_source = 400.0;  // Effective temperature of the radiation source
    double lambda = 1e-6;  // Wavelength in meters
    double omega_x = cos(30 * M_PI / 180.0);
    double omega_y = sin(30 * M_PI / 180.0);





    for (int i = 0; i < 3; i++) {
        for (int j = APERTURE_START; j < APERTURE_END; j++) {
            // Compute Planck source intensity for 6000K at 30 degrees
            double source_term = (2.0 * h * pow(c, 2)) / pow(lambda, 5) * (1.0 / (exp((h * c) / (lambda * k_B * T_source)) - 1.0));
            double weight_sum = 0.0;
            double weighted_intensity = 0.0;
            // Compute weighted sum over all discrete angles
            for (int n = 0; n < ANGLES; n++) {
                double weight = fabs(cos((theta[n] - 30) * M_PI / 180.0));  // Weight based on closeness to 30°
                weighted_intensity = weight * source_term;
                printf("I i j weight %d %d %f %f %f\n", i, j, weight, weighted_intensity,weight*planck_radiance(T_source));
                 I[i][j][n] += weighted_intensity;
                
            }

            /*for (int n = 0; n < ANGLES; n++) {
                double weight = fabs(cos((theta[n] - 30) * M_PI / 180.0));  // Weight based on closeness to 30°
                weighted_intensity += weight * source_term/weight_sum;
                if(weighted_intensity > 0) {
                    I[i][j][n] += weighted_intensity;
                } else {
                    I[i][j][n] = 0.0;  // No radiation case
                }
                
            }*/


            
        }
    }
}

void compute_average_intensity(double I[NX][NY][ANGLES], double I_avg[NX][NY]) {
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            I_avg[i][j] = 0.0; // Initialize average intensity

            // Sum over all discrete angles
            for (int n = 0; n < ANGLES; n++) {
                I_avg[i][j] += I[i][j][n];
            }

            // Compute average intensity
            I_avg[i][j] /= ANGLES;
        }
    }
}

void apply_reflective_boundaries(double I[NX][NY][ANGLES]) {
    for (int i = 0; i < NX; i++) {
        for (int n = 0; n < ANGLES; n++) {
            double omega_y = sin(theta[n] * M_PI / 180.0);

            // Reflect left boundary
            if (omega_y < 0) {  // Incoming radiation
                I[i][0][n] = I[i][1][n];  // Mirror intensity
            }

            // Reflect right boundary
            if (omega_y > 0) {  // Incoming radiation
                I[i][NY-1][n] = I[i][NY-2][n];  // Mirror intensity
            }
        }
    }
}


// Multi-angle Crank-Nicholson solver using Hypre
void solve_multi_angle_crank_nicholson(double T[NX][NY]) {
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    HYPRE_Solver solver;

    int ilower[2] = {0, 0};
    int iupper[2] = {NX-1, NY-1};

    // Initialize Hypre grid
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);

    // Create stencil for SN transport
    int stencil_size = 5;
    HYPRE_StructStencilCreate(2, stencil_size, &stencil);
    int offsets[5][2] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    for (int i = 0; i < stencil_size; i++)
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);

    // Create system matrix A
    //HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, &A);
    //HYPRE_StructMatrixSetStencil(A, stencil);
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    // Create vectors for solving the implicit system
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);


//introduce_source_term(I);
    // Time loop for Crank-Nicholson solving
        for (int t = 0; t < TIME_STEPS; t++) {
            printf("Solving , time step %d...\n", t);
    // Loop over multiple discrete angles
    for (int n = 0; n < ANGLES; n++) {
         printf("Solving , angle %d...\n", n);
        double omega_x = cos(theta[n] * M_PI / 180.0);
        double omega_y = sin(theta[n] * M_PI / 180.0);

        // Define multi-angle transport matrix coefficients
        double values[5] = {1.0 + DT * ALPHA / 2, -ALPHA / 2 * omega_x, -ALPHA / 2 * omega_x, -ALPHA / 2 * omega_y, -ALPHA / 2 * omega_y};
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++) {
                int index[2] = {i, j};
                HYPRE_StructMatrixSetValues(A, index, stencil_size, offsets,values);
            }
        HYPRE_StructMatrixAssemble(A);
        
    
            
            // Compute RHS vector (source term)
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int index[2] = {i, j};
                    double source_term = planck_radiance(T[i][j]);

                    double rhs_value = 0.0001* I[i][j][n] - DT * ALPHA / 2 * I[i][j][n] + DT * source_term;
                    HYPRE_StructVectorSetValues(b, index, rhs_value);
                }
            }
            HYPRE_StructVectorAssemble(b);

            // Solve using GMRES Krylov solver
            HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver);
            HYPRE_StructGMRESSetTol(solver, 1e-6);
            HYPRE_StructGMRESSetMaxIter(solver, 100);
            HYPRE_StructGMRESSetup(solver, A, b, x);
            HYPRE_StructGMRESSolve(solver, A, b, x);

            // Update radiation field per angle using solved values
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int index[2] = {i, j};
                    double sol_value;
                    HYPRE_StructVectorGetValues(x, index, &sol_value);
                    I[i][j][n] = sol_value;  // Store intensity for current angle
                }
            }

            // Destroy solver for reuse in next timestep
            HYPRE_StructGMRESDestroy(solver);
        }

    apply_reflective_boundaries(I);

    // Compute average intensity across all angles
    double J[NX][NY];
    compute_average_intensity(I, J) ;

    // Update temperature field using average intensity
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++) {
            T[i][j] += J[i][j]*sigmaa(T[i][j])*heat_capacity(T[i][j]) * DT;  // Update temperature based on average intensity
        }
    if(t%10==0)
            write_vtk_file(T, t/10);

    }




    

    // Cleanup Hypre memory
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
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




// Main function
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    // Initialize temperature field (room temperature)
    double T[NX][NY];
    for (int i = 0; i < NX; i++)
        for (int j = 0; j < NY; j++)
            T[i][j] = 300.0;

    // Solve transport using Crank-Nicholson method for multiple angles
    solve_multi_angle_crank_nicholson(T);

    MPI_Finalize();
    return 0;
}