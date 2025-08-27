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
#include "include/setup.h"
#include "include/boundary.h"
#include "include/material.h"
#include "include/physics.h"
#include "include/solver.h"
#include "include/utilities.h"

//Using the GMRES FE solver set the system as a top hat configuration
// 2 rows of 5 square boxes
// mats v. high heat capcity, v high aborption, low emissivity
//      v low heat capacity  v moderate absorption, higher emission

//set in top hot configuration

//the diffusion model is set in the "noddy method"

//corrected vtk output writer







int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    HYPRE_Init();

    bool dtconverged=false;
    double t=0.0;
    const int n = NX * NY * NZ;
    std::vector<double> E(n), E_new(n), source(n), cv(n), mu_a(n);
    double Ea[NX][NY][NZ]; // 3D array for temperature
    double tnp1mtn[NX][NY][NZ]; // 3D array for temperature difference
    double Tc[NX][NY][NZ]; // 3D array for temperature at current time step
    double Tn[NX][NY][NZ]; // 3D array for temperature at current time step
    double sum1, sum2, sum3,sum4;
    double Eag[num_freq_bins][NX][NY][NZ]; // 4D array for temperature
    double GradEg[num_freq_bins][NX][NY][NZ][3]={}; // gradient of the energy
    double diffusion_coeff[num_freq_bins][NX][NY][NZ]; // corrected diffusion coefficients
    double ddelr[num_freq_bins][NX][NY][NZ]; // the divergence of the (energy gradient multiplied by the corrected diffusion coefficient)
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
            GradEg[n][i][j][k][3]=0; // gradient of the energy
            diffusion_coeff[n][i][j][k]=0; // corrected diffusion coefficients
            ddelr[n][i][j][k]=0; // the divergence of the (energy gradient multiplied by the corrected diffusion coefficient)

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
    int nconverged=0;  //number of steps for which converged 
    for (int step = 0; step < NSTEP; ++step) {

        int niter=0;
        
        while(!dtconverged && step<NSTEP) {
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
        gradenergy(GradEg, Eag); // Compute the gradient of the energy
        
        //write_vtk_file(Bagn[0], 0);

        std::vector<double> values(7);

        double a_nu = 0.7; // Initialize a_nu  should be computed using rosseland mean opacity
        double kappa_nu ;
        double sumrhs=0.0;
        double sumcdt=0.0;
        double sumd=0.0;
        double summu=0.0;
        double sumdiag=0.0;
        double sumemis=0.0;
        for(int nf=0; nf<num_freq_bins; nf++) {
            int n= shuffled[nf] - 1; // Get the shuffled frequency index (0-based)
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {

        int idx = index(i, j, k);
        sum1=0;
        sum2=cv[index(i,j,k)];
        sum3=0.0;
        sum4=0.0;
        for(int nf1=0; nf1<num_freq_bins; nf1++) {
        
            sum1+=c*mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(nf1+1)*1.0e14);
            sum3=c*mu_a[index(i,j,k)]*(Bag[nf1][i][j][k]);
            sum4=sum4+c*mu_a[index(i,j,k)]*Eag[nf1][i][j][k];
            diffusion_coeff[n][i][j][k]=larsendelimiter(mu_a[idx], Eag,GradEg,i,j,k,nf1);
            ddelr[n][i][j][k]=divergence(diffusion_coeff,GradEg,i,j,k,nf1);
        }
        kappa_nu=1.0/(sum2+sum1);

            
            //double D = 1.0 / (3.0 * mu_a[idx]);
            
            //double diag = SCALE*((1.0/ (c*dt)) + 6.0 * D / (dx * dy) + mu_a[idx]);
            double diag = SCALE*((1.0/ (c*dt)) + ddelr[n][i][j][k] + mu_a[idx]);
            double D=ddelr[n][i][j][k];
            values[0] = diag;
            for (int s = 1; s < 7; ++s) values[s] = -D ;

            int ijk[3] = {i, j, k};
            // old HYPRE_StructMatrixSetValues(A, ijk, 7, (int[]){0,1,2,3,4,5,6}, values.data());
             HYPRE_StructMatrixSetValues(A, ijk, 7, (HYPRE_Int*)iindex, values.data());
             
            //double emission = mu_a[index(i,j,k)]*dBnudT(Tc[i][j][k],(nf+1)*1.0e14)*kappa_nu*dt*(sum4    -sum3);
            //double emission = mu_a[index(i,j,k)] * (Bag[n][i][j][k] - Eag[n][i][j][k]); // classic source term
            double emission = EMISSCALE*mu_a[index(i,j,k)] * Bag[n][i][j][k]; // if treating B_nu as external source
            //emission=0.0;
            double rhs = (Eag[n][i][j][k]  / (c*dt)) + mu_a[idx]*Bag[n][i][j][k] +emission;

            sumrhs+=rhs;
            sumemis=sumemis+= emission;
            sumdiag+=diag;
            sumcdt+=SCALE*(1.0/ (c*dt));
            summu+=SCALE*mu_a[idx];
            sumd+=SCALE*6.0 * D / (dx * dx);
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
        double delTmax=0.0;
        double delT;
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
                delT=fabs(Tn[i][j][k]-Tc[i][j][k]);
                if(delT>delTmax) delTmax=delT;
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
                    refg=0.0;
                    int idx = index(i, j, k);

                    for(int n=0; n<num_freq_bins; n++) 
                    {
                        emission= B_nu(Tn[i][j][k], (n+1)*1.0e14); // Example frequency bin
                        //D=1.0 / (3.0 * mu_a[idx]);
                        D=diffusion_coeff[n][i][j][k];
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
                           Tn[i][j][k]=1.0e8*std::exp(-((i-2)*(i-2)+(j-(NY)/2)*(j-(NY)/2))/4)*TINI; // top hot configuration
                }
            }
                niter++;


                if(delTmax<temptol   || niter>maxtiter)
                {
                    for (int k = 0; k < NZ; ++k)
                    for (int j = 0; j < NY; ++j)
                    for (int i = 0; i < NX; ++i) {
                        dtconverged=true;               
                        Tc[i][j][k] = std::max(Tn[i][j][k], TMIN); // Ensure non-negative temperature
                        for(int n=0; n<num_freq_bins; n++) 
                        {
                            Bag[n][i][j][k] = B_nu(Tc[i][j][k], (n+1)*1.0e14); // Update Bag for next iteration
                            Eag[n][i][j][k] = Eagn[n][i][j][k]; // Update Eag for next iteration
                        }
                    }

                    nconverged++;

                }
                else
                {
                    if(dt>dtmin)
                        dt=dt*0.7;
                    nconverged=0;

                }

                if(nconverged>3)
                {
                    if(dt<dtmax)    
                        dt=dt*1.3;
                    nconverged=0;
                }

               
               if(step%50==0)
               {
                if(dtconverged)
               {
                    std::cout << "DT CONVERGED " << step <<"   " << nconverged<<"  " <<niter <<"   "  << dt<<"  "   <<  t  << std::endl;
                    t+=dt;
                }
                else
                    std::cout << "DT NOT CONVERGED " << step << "  "<< nconverged << "   "  <<niter <<"   "  << dt<<"  "   <<  t << std::endl; 
               }
        


        //manage conservation of energy and shift energy
        double deltae;
        if(dtconverged)
        {
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
            }
  
        //compute final temperature
        // Write VTK output for visualization
        if(dtconverged && step%N_SAVEINTERVAL==0)
            write_vtk_file(Tc, (step/N_SAVEINTERVAL));

        if(step%N_SAVEINTERVAL==0)
        {
            if(dtconverged)
            {
                std::cout << "DT CONVERGED " <<"   "  <<niter <<"   "  << dt<<"  "   <<  t  << std::endl;
                t+=dt;
            }
            else
                std::cout << "DT NOT CONVERGED " << std::endl;  
            std::cout << "Step: " << step << ", Total Energy: " << etotn <<  "   " <<etot<< "   " << etotn-etot <<std::endl;
            std::cout << "Sum of RHS: " << sumrhs/(NX*NY*NZ*num_freq_bins) << ", Sum of Diagonal: " << sumdiag/(NX*NY*NZ*num_freq_bins) << ", Sum of Emission: " << sumemis/(NX*NY*NZ*num_freq_bins) << std::endl;
            std::cout << "Sum  cdt: " << sumcdt/(NX*NY*NZ*num_freq_bins) << ", Sum mu: " << summu/(NX*NY*NZ*num_freq_bins) << ", Sum D: " << sumd/(NX*NY*NZ*num_freq_bins) << std::endl;
            std::cout << "Average Temperature: " << etotn/(NX*NY*NZ*num_freq_bins) << std::endl;
        }



    }//dtnotcoverged loop
        dtconverged=false; 
        
        
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