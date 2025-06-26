#include <HYPRE_struct_ls.h>
#include <mpi.h>
#include <cmath>
#include <vector>
#include "ex.h"

//original version

constexpr int NX = 64, NY = 160, NZ = 1;
constexpr int NSTEP = 100;
constexpr double dx = 1.0 / NX;
constexpr double dt = 0.01;  // time step size

double absorption_coeff(int i) {
    return 0.1 + 10.0 * (double(i) / NX);
}

double specific_heat(int j) {
    return 1.0 + 2.0 * (double(j) / NY);
}

double planck_source(double T) {
    return 4.0 * 5.67e-8 * std::pow(T, 4);  // Simplified gray approx.
}

double initial_temperature(int i) {
    return (i == 0) ? 1000.0 : 300.0;
}

int index(int i, int j, int k) {
    return i + NX * (j + NY * k);
}


void write_vtk_file(double T[NX][NY][NZ], int timestep) {
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
    fprintf(vtkFile, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
    fprintf(vtkFile, "ORIGIN 0 0 0\n");
    fprintf(vtkFile, "SPACING 1 1 1\n");

    // Write data section
    fprintf(vtkFile, "POINT_DATA %d\n", (NX * NY*NZ));
    fprintf(vtkFile, "SCALARS Temperature float 1\n");
    fprintf(vtkFile, "LOOKUP_TABLE default\n");
 
    // Write temperature values in row-major order
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++)
        for (int k = 0; k < NZ; k++) {
            fprintf(vtkFile, "%f ", T[i][j][k]);
        }
        fprintf(vtkFile, "\n");
    }

    fclose(vtkFile);
    printf("VTK file written successfully: %s\n", filename);
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    HYPRE_Init();

    const int n = NX * NY * NZ;
    std::vector<double> E(n), E_new(n), source(n), cv(n), mu_a(n);
    double Ea[NX][NY][NZ]; // 3D array for temperature
    // Initialize material properties and initial E
    for (int k = 0; k < NZ; ++k)
    for (int j = 0; j < NY; ++j)
    for (int i = 0; i < NX; ++i) {
        int idx = index(i, j, k);
        double T = initial_temperature(i);
        E[idx] = planck_source(T);
        mu_a[idx] = absorption_coeff(i);
        cv[idx] = specific_heat(j);
        source[idx] = planck_source(T);
    }

    // HYPRE grid, stencil, matrix, vectors
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    HYPRE_StructSolver solver;

    // --- Grid definition ---
    HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);
    int ilower[3] = {0, 0, 0}, iupper[3] = {NX - 1, NY - 1, NZ - 1};
    HYPRE_StructGridSetExtents(grid, ilower, iupper);
    HYPRE_StructGridAssemble(grid);

    // --- Stencil definition ---
    HYPRE_StructStencilCreate(3, 7, &stencil);
    int offsets[7][3] = {{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    for (int s = 0; s < 7; ++s)
        HYPRE_StructStencilSetElement(stencil, s, offsets[s]);

    // --- Create and initialize matrix and vectors ---
    HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
    HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);
    
    int iindex[7]={0,1,2,3,4,5,6};
    // --- Time loop ---
    for (int step = 0; step < NSTEP; ++step) {
        std::vector<double> values(7);
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {
            int idx = index(i, j, k);
            double D = 1.0 / (3.0 * mu_a[idx]);
            double diag = cv[idx] / dt + 6.0 * D / (dx * dx) + mu_a[idx];

            values[0] = diag;
            for (int s = 1; s < 7; ++s) values[s] = -D / (dx * dx);

            int ijk[3] = {i, j, k};
            // old HYPRE_StructMatrixSetValues(A, ijk, 7, (int[]){0,1,2,3,4,5,6}, values.data());
             HYPRE_StructMatrixSetValues(A, ijk, 7, (HYPRE_Int*)iindex, values.data());
             

            double rhs = E[idx] * cv[idx] / dt + source[idx];
            //HYPRE_StructVectorSetValues(b, ijk, &rhs);
            HYPRE_StructVectorSetValues(b, ijk, rhs);
        }

        HYPRE_StructMatrixAssemble(A);
        HYPRE_StructVectorAssemble(b);
        HYPRE_StructVectorAssemble(x);

        HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructSMGSetup(solver, A, b, x);
        HYPRE_StructSMGSolve(solver, A, b, x);
        HYPRE_StructSMGDestroy(solver);

        // Extract solution
        HYPRE_StructVectorGetBoxValues(x, ilower, iupper, E_new.data());
        E = E_new;
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {
            Ea[i][j][k] = E[index(i, j, k)];
        }

                // Write VTK output for visualization
    
        if(step%10==0)
            write_vtk_file(Ea, step/10);




    }

    // Cleanup
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructGridDestroy(grid);
    HYPRE_Finalize();
    MPI_Finalize();
    return 0;
}


/*
//Krylov solver plugin for HYPRE

// Create Krylov solver and preconditioner
HYPRE_StructSolver krylov_solver, precond;
HYPRE_StructPCGCreate(MPI_COMM_WORLD, &krylov_solver);
HYPRE_StructDiagScaleCreate(MPI_COMM_WORLD, &precond);

// Set preconditioner and parameters
HYPRE_StructPCGSetPrecond(krylov_solver,
                          (HYPRE_PtrToStructSolverFcn) HYPRE_StructDiagScaleSolve,
                          (HYPRE_PtrToStructSolverFcn) HYPRE_StructDiagScaleSetup,
                          precond);

HYPRE_StructPCGSetMaxIter(krylov_solver, 200);
HYPRE_StructPCGSetTol(krylov_solver, 1e-8);
HYPRE_StructPCGSetLogging(krylov_solver, 1);

// Set up and solve
HYPRE_StructPCGSetup(krylov_solver, A, b, x);
HYPRE_StructPCGSolve(krylov_solver, A, b, x);

// Clean up
HYPRE_StructDiagScaleDestroy(precond);
HYPRE_StructPCGDestroy(krylov_solver);


*/