#include "HYPRE.h"
#include "HYPRE_krylov.h"
#include <HYPRE_struct_ls.h>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <algorithm> // for std::shuffle
#include <random> // for std::mt19937 and std::random_device
#include "ex.h"


//This version is developed from version 6
//but we attempt a more rigorous approach to the temperature iteration
//with linearization over multiple frequency groups

//Using the GMRES FE solver set the system as a top hat configuration
// 2 rows of 5 square boxes
// mats v. high heat capcity, v high aborption, low emissivity
//      v low heat capacity  v moderate absorption, higher emission

//set in top hot configuration

//the diffusion model is set in the "noddy method"

//corrected vtk output writer




constexpr int NX = 160, NY = 64, NZ = 1;
constexpr int NSTEP = 1000;
constexpr double dx = 1.0 / NX;
constexpr double dt = 0.1;  // time step size
const int num_freq_bins = 10; // Number of frequency bins
const double c = 3e10; // Speed of light in cm/s
const double a = 7.5646e-15; // Radiation constant in erg/cm^3/K^4
const double sigma = 5.6704e-5; // Stefan-Boltzmann constant in erg/cm^2/s/K^4
const double h = 6.626e-27; // Planck's constant in erg.s
const double k = 1.38064852e-16; // Boltzmann constant in erg/K
const double TMIN = 293.0; // Minimum temperature - comfortable room temperature
const double TMAX = 10000.0; // Maximum temperature

const double temptol=0.1;
const int maxtiter=10;  //temperature iteration

int boundary(int i, int j) {
    // Define the boundaries based on the grid indices
    // This function returns an integer representing the boundary type:
    int ibound=0;
    if (i < NX / 5 && (j == NY / 2   || j==0)) {
        //set up and lower wall to initial condition
        ibound = 1; // left box
    } else if (i <= (NX / 5) && j > NY / 2) {
        //left hand wall
        ibound = 2; // right box
    } else if (i > ( NX / 5) && i<(2*NX/5) && j ==0) {
        //lower boundary 2nd box
        ibound = 1; // right box
    } else if (i > (NX / 5) && i < (2 * NX / 5) && j == (NY -1)) {
        //2nd box upper section
        ibound = 2; // middle box
    }  else if (i >= (2*NX / 5) && i < (3 * NX / 5) && j < (NY/2 )) {
        //middle box lower section
        ibound = 2; // middle box
    } else if (i >= (2*NX / 5) && i < (3 * NX / 5) && j == (NY-1 )) {
        //middle box upper section
        ibound = 2; // middle box
    }  else if (i >= (3*NX / 5) && i < (4 * NX / 5) && j == 0) {
        // box 4 low section
        ibound = 2; // middle box
    }  else if (i >= (3*NX / 5) && i < (4 * NX / 5) && j == (NY-1)) {
        // box 4 low section
        ibound = 2; // middle box
    }  else if (i >= (4*NX / 5) && i < ( NX ) && j == 0) {
        // box 4 low section
        ibound = 2; // middle box
    }  else if (i >= (4*NX / 5) && i < ( NX ) && j > (NY/2)) {
        // box 5 upper section
        ibound = 2; // middle box
    }    
    else {
        ibound = 0; // background
    }


    // Return the boundary type
    // 0 = background, 1 = left box, 2 = right box,
    return ibound;
}


double absorption_coeff(int i, int j, int freq_bin = 0) {
    //frequency dependent absorption coefficient
    //for this example we use a constant value for simplicity
    double siga=1.0;
    double sigaupper=0.00001;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        siga = sigaupper ;
    } else {
        siga = 0.5;
    }
    
    return siga;
}

double specific_heat(int i, int j, int freq_bin = 0) {
    //frequency dependent heat capacity
    //for this example we use a constant value for simplicity
    double ca=200.0;
    double caupper=20000.0;

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
    if(i<4 &&  j< 3*(NY/2)/4   && j>NY/4) {
    //if(i<4 &&  j< (NY/2)/2   ) {
        tini=10000.0; // top hot configuration
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
    double sum1, sum2, sum3;
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
        if (boundary(i, j) == 1 || boundary(i, j) == 2 || boundary(i, j) == 3) {
                    T = initial_temperature(i,j); // Reset temperature for boundary conditions
                }
        Tc[i][j][k] = T; // Initialize temperature at current time step
        Tn[i][j][k] = T; // new temp
        tnp1mtn[i][j][k] = 0.0; // Initialize temperature difference
        E[idx] = planck_source(T);

        source[idx] = planck_source(T);
        for(int n=0; n<num_freq_bins; n++) {
            mu_a[idx] = absorption_coeff(i,j,n);
            cv[idx] = specific_heat(i,j,n);
            Eag[n][i][j][k] = B_nu(1.01*T, (n+1)*1e14); // Example frequency bin
            Bag[n][i][j][k] = B_nu(T, (n+1)*1e14); // Example frequency bin
            Eagn[n][i][j][k] = B_nu(1.01*T, (n+1)*1e14); // Example frequency bin
            Bagn[n][i][j][k] = B_nu(T, (n+1)*1e14); // Example frequency bin

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
            int idx = index(i, j, k);
            Eagn[0][i][j][k] = cv[idx];
            
            //std::cout << i << "  " << j  << "   " <<  k << "  " <<  index(i,j,k) <<   "  " << Ea[i][j][k] << std::endl;
        }
    write_vtk_file(Eagn[0], 0);


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

    std::vector<int> ordered;
    for (int i = 1; i <= num_freq_bins; ++i) {
    ordered.push_back(i);
    }

    // Create a random engine seeded with a non-deterministic random device
    std::random_device rd;
    std::mt19937 g(rd());


    // --- Time loop ---
    for (int step = 0; step < NSTEP; ++step) {


        // Shuffle the vector of indices for the frequencies // this reduces numerical noise
        std::vector<int> shuffled = ordered;
        std::shuffle(shuffled.begin(), shuffled.end(), g);

        // Print the shuffled vector
        /*for (int n : shuffled) {
        std::cout << n << " ";
        }
        std::cout << std::endl;*/
        //start iteration loop

        //compute tnp1-tn
        
            for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {
                int idx = index(i, j, k);
                sum1=0;
                sum2=cv[index(i,j,k)]/dt;
                sum3=0.0;
                for(int n=0; n<num_freq_bins; n++) {
                    mu_a[idx] = absorption_coeff(i,j,n);
                    cv[idx] = specific_heat(i,j,n);
                    Bag[n][i][j][k]=Bagn[n][i][j][k];                
                    sum1+=c*mu_a[index(i,j,k)]*(Eagn[n][i][j][k]-Bag[n][i][j][k]);
                    sum2+=c*mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(n+1)*1.0e14);
                    sum3+=B_nu(Tc[i][j][k],(n+1)*1.0e14) * mu_a[index(i,j,k)];

            }
            source[idx] = sum3; // Update source term
            tnp1mtn[i][j][k]=sum1/sum2;
        }



        std::vector<double> values(7);


        for(int nf=0; nf<num_freq_bins; nf++) {
            int n= shuffled[nf] - 1; // Get the shuffled frequency index (0-based)
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
            //Ea[i][j][k] = E[index(i, j, k)];
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
                //apply boundary condition


                if (boundary(i, j) == 1 || boundary(i, j) == 2 || boundary(i, j) == 3) {
                    Tn[i][j][k] = initial_temperature(i,j); // Reset temperature for boundary conditions
                }

                Tc[i][j][k] = std::max(Tn[i][j][k], TMIN); // Ensure non-negative temperature


                for(int n=0; n<num_freq_bins; n++) 
                {
                    Bagn[n][i][j][k] = B_nu(Tc[i][j][k], (n+1)*1.0e14); // Update Bag for next iteration
                }
        }

        //compute final temperature
        // Write VTK output for visualization
        if(step%100==0)
            write_vtk_file(Tc, 1+(step/100));




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