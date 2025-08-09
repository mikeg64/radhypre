#include "HYPRE.h"
#include "HYPRE_krylov.h"
#include <HYPRE_struct_ls.h>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <algorithm> // for std::shuffle
#include <random> // for std::mt19937 and std::random_device
#include <functional>
#include <iostream>

#include "ex.h"

//Using the GMRES FE solver set the system as a top hat configuration
// 2 rows of 5 square boxes
// mats v. high heat capcity, v high aborption, low emissivity
//      v low heat capacity  v moderate absorption, higher emission

//set in top hot configuration

//the diffusion model is set in the "noddy method"

//corrected vtk output writer




constexpr int NX = 160, NY = 64, NZ = 1;
constexpr int NSTEP = 100000;
constexpr int N_SAVEINTERVAL=100;
constexpr double dx = 0.07 / NX;
constexpr double dy = 0.05 / NY;
constexpr double dt = 0.01;  // time step size
const int num_freq_bins = 10; // Number of frequency bins
const double c = 3e10; // Speed of light in cm/s
const double a = 7.5646e-15; // Radiation constant in erg/cm^3/K^4
const double sigma = 5.6704e-5; // Stefan-Boltzmann constant in erg/cm^2/s/K^4
const double h = 6.626e-27; // Planck's constant in erg.s
const double k = 1.38064852e-16; // Boltzmann constant in erg/K
const double TMIN = 293.0; // Minimum temperature - comfortable room temperature
const double TMAX = 10000.0; // Maximum temperature
const double TINI = 1000000.0;//11604525.0061598; // Maximum temperature
const double EINILT2 = 0.000000000001; // Initial temperature for top hot configuration
const double EINIEQ2 = 0.0000000001; // Initial temperature for equilibrium configuration
const int    BOUNDTYPE = 2; // Boundary type for the fixed temp is 1 reflected energy is 2, and background is 0, 3 do nothing
const double temptol=0.1;
const int maxtiter=10;  //temperature iteration
const int BUY=1;//upper y 1
const int BLY=2;//lower y 2
const int BLX=3;//left x 3
const int BRX=4;//right x 4

int boundary(int i, int j) {
    // Define the boundaries based on the grid indices
    // This function returns an integer representing the boundary type:
    int ibound=0;
    if (i < NX / 5 && (j == NY / 2   || j==0)) {
        //set up and lower wall to initial condition
        ibound = BOUNDTYPE; // left box
    } else if (i <= (NX / 5) && j > NY / 2) {
        //left hand wall
        ibound = BOUNDTYPE; // right box
    } else if (i > ( NX / 5) && i<(2*NX/5) && j ==0) {
        //lower boundary 2nd box
        ibound = BOUNDTYPE; // right box
    } else if (i > (NX / 5) && i < (2 * NX / 5) && j == (NY -1)) {
        //2nd box upper section
        ibound = BOUNDTYPE; // middle box
    }  else if (i >= (2*NX / 5) && i < (3 * NX / 5) && j < (NY/2 )) {
        //middle box lower section
        ibound = BOUNDTYPE; // middle box
    } else if (i >= (2*NX / 5) && i < (3 * NX / 5) && j == (NY-1 )) {
        //middle box upper section
        ibound = BOUNDTYPE; // middle box
    }  else if (i >= (3*NX / 5) && i < (4 * NX / 5) && j == 0) {
        // box 4 low section
        ibound = BOUNDTYPE; // middle box
    }  else if (i >= (3*NX / 5) && i < (4 * NX / 5) && j == (NY-1)) {
        // box 4 low section
        ibound = BOUNDTYPE; // middle box
    }  else if (i >= (4*NX / 5) && i < ( NX ) && j == 0) {
        // box 4 low section
        ibound = BOUNDTYPE; // middle box
    }  else if (i >= (4*NX / 5) && i < ( NX ) && j > (NY/2)) {
        // box 5 upper section
        ibound = BOUNDTYPE; // middle box
    }    
    else {
        ibound = 0; // background
    }


    // Return the boundary type
    // 0 = background, 1 = left box, 2 = right box,
    return ibound;
}


//used for the reflected energy boundary condition
//upper y 1  BUY
//lower y 2   BLY
//left x 3     BLX
//right x 4    BRX
int refboundarytype(int i, int j) {
    // Define the boundaries based on the grid indices
    // This function returns an integer representing the boundary type:
    int ibound=0;
    if (i < NX / 5 && (j == NY / 2   || j==0)) {
        //set up and lower wall to initial condition
        ibound = (j==0)*BLY+(j == NY / 2)*BUY; // left box
    } else if (i <= (NX / 5) && j > NY / 2) {
        //left hand wall
        ibound = BLX; // right box
    } else if (i > ( NX / 5) && i<(2*NX/5) && j ==0) {
        //lower boundary 2nd box
        ibound = (BLY); // right box
    } else if (i > (NX / 5) && i < (2 * NX / 5) && j == (NY -1)) {
        //2nd box upper section
        ibound = BUY; // middle box
    }  else if (i >= (2*NX / 5) && i < (3 * NX / 5) && j < (NY/2 )) {
        //middle box lower section
        ibound = BLY; // middle box
    } else if (i >= (2*NX / 5) && i < (3 * NX / 5) && j == (NY-1 )) {
        //middle box upper section
        ibound = BUY; // middle box
    }  else if (i >= (3*NX / 5) && i < (4 * NX / 5) && j == 0) {
        // box 4 low section
        ibound = BLY; // middle box
    }  else if (i >= (3*NX / 5) && i < (4 * NX / 5) && j == (NY-1)) {
        // box 4 low section
        ibound = BUY; // middle box
    }  else if (i >= (4*NX / 5) && i < ( NX ) && j==0) {
        // box 4 low section
        ibound = BLY; // middle box
    }  else if (i >= (4*NX / 5) && i < ( NX ) && j > (NY/2)) {
        // box 5 upper section
        ibound = BUY; // middle box
    }    
    else {
        ibound = 0; // background
    }


    // Return the boundary type
    // 0 = background, 1 = left box, 2 = right box,
    return ibound;
}

double ereflect(int n, int i, int j)
{
    double eref=0.0;

    return eref;
}

double absorption_coeff(int i, int j) {
    double siga;
    double sigaupper=200.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        siga = sigaupper ;
    } else {
        siga = 0.5;
    }
    
    return siga;
}

double specific_heat(int i, int j) {
    double ca;
    double caupper=10000.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        ca = caupper ;
    } else {
        ca = 1;
    }
    
    return ca;
}

//same as B_nu but for Planck's law
// This function computes the spectral radiance of a black body at frequency nu and temperature T
double planck_emission(double nu, double T) {
    //const double h = 6.62607015e-27;     // erg·s
    //const double c = 2.99792458e10;      // cm/s
    //const double k = 1.380649e-16;       // erg/K

    double numerator = 2.0 * h * std::pow(nu, 3) / (c * c);
    double exponent = h * nu / (k * T);
    double denominator = std::exp(exponent) - 1.0;

    return numerator / denominator; // erg·s⁻¹·cm⁻²·Hz⁻¹·sr⁻¹
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


double dBdT(double nu, double T) {
    const double h = 6.62607015e-27;     // erg·s
    const double c = 2.99792458e10;      // cm/s
    const double k = 1.380649e-16;       // erg/K

    double x = h * nu / (k * T);
    double exp_x = std::exp(x);
    double factor = (2.0 * h * h * nu * nu * nu * nu) / (c * c * k * T * T);
    return factor * exp_x / std::pow(exp_x - 1.0, 2);
}


// Define the opacity function
std::function<double(double)> make_kappa_nu(double kappa0, double nu0, double alpha) {
    return [=](double nu) {
        return kappa0 * std::pow(nu / nu0, -alpha);
    };
}
// Function to compute the Rosseland mean opacity for a group of frequencies
// kappa_nu is a function that returns the opacity at frequency nu

//example usage of function
//    double T = 6000.0;                  // Temperature in Kelvin
//    double nu_min = 1e13;              // Minimum frequency in Hz
//    double nu_max = 1e15;              // Maximum frequency in Hz
//    double kappa0 = 1.0;               // Reference opacity
//    double nu0 = 1e14;                 // Reference frequency
//    double alpha = 1.5;                // Power-law index

    // Create the opacity function
//    auto kappa_nu = make_kappa_nu(kappa0, nu0, alpha);

    // Compute Rosseland mean opacity for the group
//    double kappa_R = rosseland_mean_opacity_group(kappa_nu, T, nu_min, nu_max);








double rosseland_mean_opacity_group(
    std::function<double(double)> kappa_nu, // frequency-dependent opacity
    double T,
    double nu_min,
    double nu_max,
    int N = 1000 // number of integration points
) {
    double dnu = (nu_max - nu_min) / N;
    double numerator = 0.0;
    double denominator = 0.0;

    for (int i = 0; i < N; ++i) {
        double nu = nu_min + i * dnu;
        double dB_dT = dBdT(nu, T);
        numerator += (1.0 / kappa_nu(nu)) * dB_dT * dnu;
        denominator += dB_dT * dnu;
    }

    return denominator / numerator; // Rosseland mean opacity for the group
}

double rosseland_mean_opacity_group(
    double kappa, // frequency-dependent opacity
    double T,
    int num_freq_bins = 1000 // number of integration points
) {
    double dnu = 1.0e14;
    double numerator = 0.0;
    double denominator = 0.0;

 


    for (int i = 0; i < num_freq_bins; ++i) {
        double nu = (i+1)*1e14;
        double dB_dT = dBdT(nu, T);
        numerator += (1.0 / kappa) * dB_dT * dnu;
        denominator += dB_dT * dnu;
    }

    return denominator / numerator; // Rosseland mn opacity for the group
}






double planck_source(double T) {
    return 4.0 * 5.67e-8 * std::pow(T, 4);  // Simplified gray approx.
}

double initial_temperature(int i, int j) {
    double tini=300;
    if(i<4 &&  j< 3*(NY/2)/4   && j>NY/4) {
    //if(i<4 &&  j< (NY/2)/2   ) {
        tini=TINI; // top hot configuration
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
    double sum1, sum2, sum3,sum4;
    double Eag[num_freq_bins][NX][NY][NZ]; // 4D array for temperature
    double Bag[num_freq_bins][NX][NY][NZ]; // 4D array for temperature
    double Eagn[num_freq_bins][NX][NY][NZ]; // new values
    double Bagn[num_freq_bins][NX][NY][NZ]; //new values
    double etot, etotn; // total energy at current and next time step
    // Initialize material properties and initial E
    etot = 0.0;
    etotn = 0.0;
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
        mu_a[idx] = absorption_coeff(i,j);
        cv[idx] = specific_heat(i,j);
        source[idx] = planck_source(T);
        etot += cv[idx]*Tc[i][j][k]; // Calculate total energy at current time step
        for(int n=0; n<num_freq_bins; n++) {
            Eag[n][i][j][k] = B_nu(T, (n+1)*1e14); // Example frequency bin
            Bag[n][i][j][k] = B_nu(T, (n+1)*1e14); // Example frequency bin
            Eagn[n][i][j][k] = B_nu(1.01*T, (n+1)*1e14); // Example frequency bin
            Bagn[n][i][j][k] = B_nu(1.01*T, (n+1)*1e14); // Example frequency bin
            etot+=Eag[n][i][j][k];
            etotn=etot;
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
            //Eagn[0][i][j][k] = cv[idx];  //just eriting the initial value
            //Eagn[0][i][j][k] = Tc[i][j][k]; // Initialize with temperature
            for(int n=0; n<num_freq_bins; n++) {
            if(i<1 &&  j< 3*(NY/2)/4   && j>NY/4  && n<2  || n>2 ) {
                            //if(i<4 &&  j< (NY/2)/2   ) {
                           ;//Eagn[n][i][j][k]=EINILT2; // top hot configuration
                }
             if(i<1 &&  j< 3*(NY/2)/4   && j>NY/4  && n==2 ) {
                            //if(i<4 &&  j< (NY/2)/2   ) {
                           ;//Eagn[n][i][j][k]=EINIEQ2; // top hot configuration
                }
            }
            
            //std::cout << i << "  " << j  << "   " <<  k << "  " <<  index(i,j,k) <<   "  " << Ea[i][j][k] << std::endl;
        }
    /*write_vtk_file(Bagn[0], 0);
    write_vtk_file(Bag[0], 1);
    write_vtk_file(Eagn[0], 2);
    write_vtk_file(Eag[0], 3);*/


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
            etot=0.0;
            for (int k = 0; k < NZ; ++k)
            for (int j = 0; j < NY; ++j)
            for (int i = 0; i < NX; ++i) {

                int idx = index(i, j, k);
                etot += cv[idx]*Tc[i][j][k]; // Calculate total energy at current time step
                sum1=0;
                sum2=cv[index(i,j,k)];
                sum4=0.0;
                for(int n=0; n<num_freq_bins; n++) {
                    ;//Bag[n][i][j][k]=Bagn[n][i][j][k];                
                    sum1+=c*mu_a[index(i,j,k)]*(Eagn[n][i][j][k]-Bag[n][i][j][k]);
                    sum2+=c*mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(n+1)*1.0e14);
                    //sum4+=B_nu(Tc[i][j][k],(n+1)*1.0e14) * mu_a[index(i,j,k)];
                    etot+=Eag[n][i][j][k];
                    }
            sum2 = 1.0/std::max(sum2, 1e-10); // Prevent division by zero
            source[idx] =0.0; sum3; // Update source term
            tnp1mtn[i][j][k]=sum1/sum2;
        }

        //write_vtk_file(Bagn[0], 0);

        std::vector<double> values(7);

        double a_nu = 0.7; // Initialize a_nu  should be computed using rosseland mean opacity
        double kappa_nu ;
        for(int nf=0; nf<num_freq_bins; nf++) {
            int n= shuffled[nf] - 1; // Get the shuffled frequency index (0-based)
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {


        sum1=0;
        sum2=cv[index(i,j,k)];
        sum3=0.0;
        sum4=0.0;
        for(int nf1=0; nf1<num_freq_bins; nf1++) {
           
            sum1+=c*mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(nf1+1)*1.0e14);
            sum3=c*mu_a[index(i,j,k)]*(Bag[nf1][i][j][k]);
            sum4=sum4+c*mu_a[index(i,j,k)]*Eag[nf1][i][j][k];
        }
        kappa_nu=1.0/(sum2+sum1);

            int idx = index(i, j, k);
            double D = 1.0 / (3.0 * mu_a[idx]);
            double diag = (1.0/ (c*dt)) + 6.0 * D / (dx * dx) + mu_a[idx];

            values[0] = diag;
            for (int s = 1; s < 7; ++s) values[s] = -D / (dx * dx);

            int ijk[3] = {i, j, k};
            // old HYPRE_StructMatrixSetValues(A, ijk, 7, (int[]){0,1,2,3,4,5,6}, values.data());
             HYPRE_StructMatrixSetValues(A, ijk, 7, (HYPRE_Int*)iindex, values.data());
             
            double emission = mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(nf+1)*1.0e14)*kappa_nu*dt*(sum4    -sum3);
            //emission=0.0;
            double rhs = (Eag[n][i][j][k]  / (c*dt)) + mu_a[idx]*Bag[n][i][j][k] +emission;
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

            /*if(i<4 &&  j< 3*(NY/2)/4   && j>NY/4  && n<2  || n>2 ){
                            //if(i<4 &&  j< (NY/2)/2   ) {

                            Eagn[n][i][j][k]= EINILT2*((1.0/(1.0-exp(-(step-6000)*dt/1.0)))-(1.0/(1.0-exp(-(step-10000)*dt/1.0)))); // top hot configuration
                }
            if(i<4 &&  j< 3*(NY/2)/4   && j>NY/4  && n==2 ) {
                            //if(i<4 &&  j< (NY/2)/2   ) {
                            Eagn[n][i][j][k]= EINIEQ2*((1.0/(1.0-exp(-(step-6000)*dt/1.0)))-(1.0/(1.0-exp(-(step-10000)*dt/1.0)))); // top hot configuration
                }  */
        }
    }  //end of loop over frequencies
        etotn=0.0;
        
        //compute new temps
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {

                
                

                sum1=0;
                sum4=cv[index(i,j,k)];
                sum2=0.0;
                for(int n=0; n<num_freq_bins; n++) 
                {
                    sum1+=c*mu_a[index(i,j,k)]*(Eagn[n][i][j][k]-Bag[n][i][j][k]);
                    sum2+=(c*mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(n+1)*1.0e14));
                    etotn+=Eagn[n][i][j][k];
                              
                }
                //Tn[i][j][k]=Tc[i][j][k] + ((sum1/sum4)/(1+sum2));
                Tn[i][j][k]=Tc[i][j][k] + dt*(sum1/(sum2+sum4));
                //apply boundary condition
                etotn += cv[index(i,j,k)]*Tn[i][j][k]; // Calculate total energy at current time step

                if (boundary(i, j) == 1 ) {
                    Tn[i][j][k] = initial_temperature(i,j); // Reset temperature for boundary conditions
                }
                if (boundary(i, j) == 2 ) {


                    //  BUY=1;//upper y 1
            //const int BLY=2;//lower y 2
            //const int BLX=3;//left x 3
            //const int BRX=4;//right x 4
                    double refg,F_in, F_out, D,emission ;
                    refg=0.7;
                    int idx = index(i, j, k);

                    for(int n=0; n<num_freq_bins; n++) 
                    {
                        emission= B_nu(Tn[i][j][k], (n+1)*1.0e14); // Example frequency bin
                        D=1.0 / (3.0 * mu_a[idx]);
                       if(refboundarytype(i,j)==BUY)
                       {
                            // Estimate incident flux from interior cell
                            F_in = -D * (Eagn[n][i][j][k] - Eagn[n][i][j-1][k]) / dy;
                            // Reflected + emitted flux
                            F_out = refg * F_in + emission;
                            Eagn[n][i][j][k] = Eagn[n][i][j][k] + F_out * dy / D;
                       }
                       else if(refboundarytype(i,j)==BLY)
                       {
                            // Estimate incident flux from interior cell
                            F_in = -D * (Eagn[n][i][j+1][k] - Eagn[n][i][j][k]) / dy;
                            // Reflected + emitted flux
                            F_out = refg * F_in + emission;
                            Eagn[n][i][j][k] = Eagn[n][i][j][k] + F_out * dy / D;
                       }
                       else if(refboundarytype(i,j)==BLX)
                       {
                             // Estimate incident flux from interior cell
                            F_in = -D * (Eagn[n][i+1][j][k] - Eagn[n][i][j][k]) / dx;
                            // Reflected + emitted flux
                            F_out = refg * F_in + emission;
                            Eagn[n][i][j][k] = Eagn[n][i][j][k] + F_out * dx / D;
                       }
                       else if(refboundarytype(i,j)==BRX)
                       {
                            // Estimate incident flux from interior cell
                            F_in = -D * (Eagn[n][i][j][k] - Eagn[n][i-1][j][k]) / dx;
                            // Reflected + emitted flux
                            F_out = refg * F_in + emission;
                            Eagn[n][i][j][k] = Eagn[n][i][j][k] + F_out * dx / D;

                       }
                       
                       
                       
                        Eagn[n][i][j][k] = ereflect(n,i,j); // Reset temperature for boundary conditions
                    }
                }
                if(i<4 &&  j< 3*(NY/2)/4   && j>NY/4) {
                            //if(i<4 &&  j< (NY/2)/2   ) {
                           Tn[i][j][k]=TINI; // top hot configuration
                }


                Tc[i][j][k] = std::max(Tn[i][j][k], TMIN); // Ensure non-negative temperature


                for(int n=0; n<num_freq_bins; n++) 
                {
                    Bag[n][i][j][k] = B_nu(Tc[i][j][k], (n+1)*1.0e14); // Update Bag for next iteration
                    Eag[n][i][j][k] = Eagn[n][i][j][k]; // Update Eag for next iteration
                }
        }


        //manage conservation of energy and shift energy
        double deltae;
        deltae=(etotn-etot)/(NX*NY*NZ*num_freq_bins);
                        for(int n=0; n<num_freq_bins; n++) 
                {
                    for (int k = 0; k < NZ; ++k)
                    for (int j = 0; j < NY; ++j)
                    for (int i = 0; i < NX; ++i) {
                        int idx = index(i, j, k);
                        if((Eagn[n][i][j][k]-deltae)>0)
                            Eagn[n][i][j][k] = Eag[n][i][j][k]-deltae;
                       
                    }
                }
  
        //compute final temperature
        // Write VTK output for visualization
        if(step%N_SAVEINTERVAL==0)
            write_vtk_file(Tc, (step/N_SAVEINTERVAL));

        if(step%N_SAVEINTERVAL==0)
            std::cout << "Step: " << step << ", Total Energy: " << etotn <<  "   " <<etot<< "   " << etotn-etot <<std::endl;




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