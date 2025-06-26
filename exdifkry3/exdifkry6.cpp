#include "HYPRE.h"
#include "HYPRE_krylov.h"
#include <HYPRE_struct_ls.h>
#include <mpi.h>
#include <cmath>
#include <vector>
#include "ex.h"

//Using the GMRES FE solver set the system as a top hot configuration
// 2 rows of 5 square boxes
// mats v. high heat capcity, v high aborption, low emissivity
//      v low heat capacity  v moderate absorption, higher emission

//set in top hot configuration

//the diffusion model is set in the "noddy method"

//corrected vtk output writer

constexpr int NX = 160, NY = 64, NZ = 1;
constexpr int NSTEP = 100;
constexpr double dx = 1.0 / NX;
constexpr double dt = 0.01;  // time step size
const int num_freq_bins = 10; // Number of frequency bins
const double c = 3e10; // Speed of light in cm/s
const double a = 7.5646e-15; // Radiation constant in erg/cm^3/K^4
const double sigma = 5.6704e-5; // Stefan-Boltzmann constant in erg/cm^2/s/K^4
const double h = 6.626e-27; // Planck's constant in erg.s
const double k = 1.38064852e-16; // Boltzmann constant in erg/K

const double temptol=0.1;
const int maxtiter=10;  //temperature iteration



double absorption_coeff(int i, int j) {
    double siga=0.1;
    double sigaupper=1000.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        siga = sigaupper ;
    } else {
        siga = 0.1;
    }
    
    return siga;
}

double specific_heat(int i, int j) {
    double ca=1.0;
    double caupper=2000.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        ca = caupper ;
    } else {
        ca = 0.1;
    }
    
    return ca;
}

// Function to compute Planck distribution
double B_nu(double T, double nu) {
    return (2 * h * pow(nu, 3) / pow(c, 2)) / (exp(h * nu / (k * T)) - 1);
}

// Function to compute Planck distribution
double dBnudT(double T, double nu) {
    //return 0.0;
    return (exp(h*nu/(k*T))*2 * pow(h,2) * pow(nu, 4) / (k*pow(T,2)*pow(c, 2))) / pow(exp(h * nu / (k * T)) - 1,2);
}

double planck_source(double T) {
    return 4.0 * 5.67e-8 * std::pow(T, 4);  // Simplified gray approx.
}

double initial_temperature(int i, int j) {
    double tini=300;
    if(i<2 &&  j< (NY/2)) {
        tini=1000.0; // top hot configuration
    } else {
        tini=300.0;
    }

    return tini;
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
     for (int j = 0; j < NY; j++)
    for (int i = 0; i < NX; i++) {
       
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
    double tnp1mtn[NX][NY][NZ]; // 3D array for temperature difference
    double Tc[NX][NY][NZ]; // 3D array for temperature at current time step
    double Tn[NX][NY][NZ]; // 3D array for temperature at current time step
    double sum1, sum2;
    double Eag[num_freq_bins][NX][NY][NZ]; // 4D array for temperature
    double Bag[num_freq_bins][NX][NY][NZ]; // 4D array for temperature
    double Eagn[num_freq_bins][NX][NY][NZ]; // new values
    double Bagn[num_freq_bins][NX][NY][NZ]; //new values
    // Initialize material properties and initial E
    for (int k = 0; k < NZ; ++k)
    for (int j = 0; j < NY; ++j)
    for (int i = 0; i < NX; ++i) {
        int idx = index(i, j, k);
        double T = initial_temperature(i,j);
        Tc[i][j][k] = T; // Initialize temperature at current time step
        Tn[i][j][k] = T; // new temp
        tnp1mtn[i][j][k] = 0.0; // Initialize temperature difference
        E[idx] = planck_source(T);
        mu_a[idx] = absorption_coeff(i,j);
        cv[idx] = specific_heat(i,j);
        source[idx] = planck_source(T);
        for(int n=0; n<num_freq_bins; n++) {
            Eag[n][i][j][k] = B_nu(T, (n+1)*1e14); // Example frequency bin
            Bag[n][i][j][k] = B_nu(T, (n+1)*1e14); // Example frequency bin
            Eagn[n][i][j][k] = B_nu(1.01*T, (n+1)*1e14); // Example frequency bin
        }
      
    }

    // HYPRE grid, stencil, matrix, vectors
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    //HYPRE_StructSolver solver;
     HYPRE_StructSolver krylov_solver;  //MKG June 2025




    for (int k = 0; k < NZ; ++k)
    for (int j = 0; j < NY; ++j)
    for (int i = 0; i < NX; ++i) {
            //Ea[i][j][k] = E[index(i, j, k)];
            Ea[i][j][k] = E[index(i, j, k)];
            //std::cout << i << "  " << j  << "   " <<  k << "  " <<  index(i,j,k) <<   "  " << Ea[i][j][k] << std::endl;
        }
    write_vtk_file(Ea, 0);


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


        //start iteration loop

        //compute tnp1-tn
        
            for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {

                sum1=0;
                sum2=cv[index(i,j,k)]/dt;
                for(int n=0; n<num_freq_bins; n++) {
                
                    sum1+=c*mu_a[index(i,j,k)]*(Eagn[n][i][j][k]-Bag[n][i][j][k]);
                    sum2+=c*mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(n+1)*1.0e14);

                


            }
            tnp1mtn[i][j][k]=sum1/sum2;
        }



        std::vector<double> values(7);


        for(int n=0; n<num_freq_bins; n++) {
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
             

            double rhs = Eag[n][i][j][k] * cv[idx] / dt + source[idx];
            //HYPRE_StructVectorSetValues(b, ijk, &rhs);
            HYPRE_StructVectorSetValues(b, ijk, rhs);
        }
       
        HYPRE_StructMatrixAssemble(A);
        HYPRE_StructVectorAssemble(b);
        HYPRE_StructVectorAssemble(x);

        // Create and initialize the HYPRE krylov solver  MKG June 2025
        // Create Krylov solver and preconditioner
        // Solve using GMRES Krylov solver
        HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &krylov_solver);
        HYPRE_StructGMRESSetTol(krylov_solver, 1e-6);
        HYPRE_StructGMRESSetMaxIter(krylov_solver, 100);
 

        // Set up and solve
        HYPRE_StructGMRESSetup(krylov_solver, A, b, x);
        HYPRE_StructGMRESSolve(krylov_solver, A, b, x);
       
        // Clean up
        //HYPRE_StructDiagScaleDestroy(precond);
        HYPRE_StructGMRESDestroy(krylov_solver);
        //HYPRE_StructGMRESDestroy(precond);

        // Extract solution
        HYPRE_StructVectorGetBoxValues(x, ilower, iupper, E_new.data());
        E = E_new;
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {
            Ea[i][j][k] = E[index(i, j, k)];
            Eagn[n][i][j][k] = E[index(i, j, k)];
        }
    }  //end of loop over frequencies

        //compute new temps
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {

                sum1=0;
                sum2=cv[index(i,j,k)]/dt;
                for(int n=0; n<num_freq_bins; n++) 
                {
                    sum1+=c*mu_a[index(i,j,k)]*(Eagn[n][i][j][k]-Bag[n][i][j][k]);
                    sum2+=c*mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(n+1)*1.0e14);
                              
                }
                Tn[i][j][k]=Tc[i][j][k] + (sum1/sum2);
                    for(int n=0; n<num_freq_bins; n++) 
                {
                    Bagn[n][i][j][k] = B_nu(Tn[i][j][k], (n+1)*1.0e14); // Update Bag for next iteration
                }
        }

        //compute final temperature
        // Write VTK output for visualization
        if(step%10==0)
            write_vtk_file(Tn, 1+(step/10));




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