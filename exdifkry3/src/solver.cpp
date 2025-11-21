#include "../include/solver.h"
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"

 


/*int setupsolver() {
    // This function sets up the solver parameters and initializes necessary variables
    // It can be expanded to include more complex setup logic if needed
    // For now, it simply returns true to indicate successful setup
    return 0;
}*/

/*void solve_radiation_groups(const Mesh& mesh, State& state) {

    for (int g = 0; g < NUM_GROUPS; ++g) {

        // Setup HYPRE solver for group g

        HYPRE_StructGrid grid;
        HYPRE_StructMatrix A;
        HYPRE_StructVector b, x;
        HYPRE_StructSolver krylov_solver;

 

        // Initialize grid, matrix, vectors (omitted for brevity)

 

        // Fill matrix A and vector b using diffusion equation

        // A[i][j] = -D_g * Laplacian + sigma_a

        // b[i] = source_term[g][i]

 

        // Solve

        HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &krylov_solver);
        HYPRE_StructGMRESSetTol(krylov_solver, 1e-6);
        HYPRE_StructGMRESSetup(krylov_solver, A, b, x);
        HYPRE_StructGMRESSolve(krylov_solver, A, b, x);

 

        // Extract solution into state.radiation_flux[g]

        // (omitted for brevity)

 

        HYPRE_StructGMRESDestroy(krylov_solver);

    }

}*/


    //RadSolve::RadSolve(int mnx, int mny) : nx(mnx), ny(mny), T(mnx, std::vector<double>(mny, 300.0)) {
    RadSolve::RadSolve(int mnx, int mny, int mnz, Pars &pars) : nx(mnx), ny(mny), nz(mnz) {
        //MPI_Init(NULL, NULL);
        std::cout<<"In RadSolve constructor"<<std::endl;
        int iindex[7]={0,1,2,3,4,5,6};
        int ilower[3] = {0, 0, 0};
        int iupper[3] = {pars.nx - 1, pars.ny - 1, pars.nz-1};
        int offsets[7][3] = {{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
        HYPRE_Init();
        // --- Grid definition ---
        HYPRE_StructGridCreate(MPI_COMM_WORLD, 3, &grid);
        HYPRE_StructGridSetExtents(grid, ilower, iupper);
        HYPRE_StructGridAssemble(grid);
        std::cout<<"In RadSolve constructor Grid"<<std::endl;
        // --- Stencil definition ---
        HYPRE_StructStencilCreate(3, 7, &stencil);

        for (int s = 0; s < 7; ++s)
            HYPRE_StructStencilSetElement(stencil, s, offsets[s]);

            
        std::cout<<"In RadSolve constructor Stencil"<<std::endl;

        // --- Create and initialize matrix and vectors ---
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
        HYPRE_StructMatrixInitialize(A);

        std::cout<<"In RadSolve constructor Matrix"<<std::endl;

        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
        HYPRE_StructVectorInitialize(b);
        HYPRE_StructVectorInitialize(x);
        
        std::cout<<"In RadSolve constructor Vectors"<<std::endl;

        // Set matrix values to zero initially
        std::vector<double> values(7, 0.0);
        std::vector<double> rhs(nx * ny * nz, 0.0);
       /* for (int k = 0; k < nz; ++k)
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i) {
                    //int ilower[3] = {i, j, k};
                    //int iupper[3] = {i, j, k};  // single point update
                    //HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 7, iindex, values.data());
                    //HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, 7, offsets, values.data());
                    int ijk[3] = {i, j, k};
                    HYPRE_StructMatrixSetValues(A, ijk, 7, (HYPRE_Int*)iindex, values.data());
                    //HYPRE_StructVectorSetValues(b, ijk, &rhs[(i + j * nx + k * nx * ny)]);
                }*/
        //HYPRE_StructMatrixAssemble(A);
        //HYPRE_StructVectorSetBoxValues(b, ilower, iupper, rhs.data());
        std::cout<<"In RadSolve constructor finished mat"<<std::endl;
        // Set vector values to zero initially
        //std::vector<double> rhs(nx * ny * nz, 0.0);
        //HYPRE_StructVectorSetBoxValues(b, ilower, iupper, rhs.data());
        // HYPRE_StructVectorSetValues(b, ijk, rhs);
        //HYPRE_StructVectorAssemble(b);

        //std::vector<double> x_init(nx * ny * nz, 0.0);
        //HYPRE_StructVectorSetBoxValues(x, ilower, iupper, x_init.data());
        //HYPRE_StructVectorAssemble(x);  
        std::cout<<"HYPRE Matrix and vectors created and initialized"<<std::endl;
        // --- Create solver ---



       //delete me
       /* HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
        HYPRE_StructStencilCreate(2, 5, &stencil);
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);*/
       //delete above
    
    
        grad_energy.resize(NUM_GROUPS);
        for(int g=0;g<NUM_GROUPS;g++) {
            grad_energy[g].resize(nx*ny, std::vector<double>(3, 0.0)); // 3 for x,y,z components
        }


        Bag.resize(NUM_GROUPS, std::vector(nx*ny*nz, 0.0));
        diff_coeff.resize(NUM_GROUPS, std::vector(nx*ny*nz, 0.0));
        ddelr.resize(NUM_GROUPS, std::vector(nx*ny*nz, 0.0));  
        
        
        for (int i = 1; i <= NUM_GROUPS; ++i) {
            ordered.push_back(i);
        }

        // Create a random engine seeded with a non-deterministic random device
        std::cout<<"In RadSolve constructor before shuffle"<<std::endl;
        



    }


    RadSolve::~RadSolve() {
        HYPRE_StructMatrixDestroy(A);
        HYPRE_StructGridDestroy(grid);
        HYPRE_StructStencilDestroy(stencil);
        HYPRE_StructVectorDestroy(b);
        HYPRE_StructVectorDestroy(x);
        //MPI_Finalize();
    }

    int RadSolve::updatestate(Pars &pars, Mesh &mesh, State &state) {   
        int status=0;
        // Update the state using radiation transport solver
        solveRadiationTransport(mesh, state, pars, pars.time);

        //careful select BCs here
        apply_milne_boundary_conditions(mesh, state, pars);
        solve_material_heating(mesh, state,pars);

        return status;
    }

    // Method to return the temperature field
    /*const std::vector<std::vector<double>>& RadSolve::getTemperatureField() const {
        return T;
    }*/

    // Example function to modify the temperature field
    /*void RadSolve::setTemperature(int i, int j, double value) {
        if (i >= 0 && i < nx && j >= 0 && j < ny) {
            T[i][j] = value;
        } else {
            std::cerr << "Error: Index out of bounds!" << std::endl;
        }
    }*/



    /*void RadSolve::readMesh(const std::string& filename) {
        std::ifstream meshFile(filename);
        std::string line;
        while (std::getline(meshFile, line)) {
            if (line.find("Nodes") != std::string::npos) {
                while (std::getline(meshFile, line)) {
                    if (line.find("EndNodes") != std::string::npos) break;
                    double x, y;
                    int index;
                    sscanf(line.c_str(), "%d %lf %lf", &index, &x, &y);
                    nodes.push_back({x, y});
                }
            }
        }
        meshFile.close();
    }*/

    /*void RadSolve::setupGrid() {
        int ilower[2] = {0, 0};
        int iupper[2] = {nx-1, ny-1};
        HYPRE_StructGridSetExtents(grid, ilower, iupper);
        HYPRE_StructGridAssemble(grid);

        int offsets[5][2] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};
        for (int i = 0; i < 5; i++)
            HYPRE_StructStencilSetElement(stencil, i, offsets[i]);

           // Create vectors
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
        HYPRE_StructVectorInitialize(b);
        HYPRE_StructVectorInitialize(x);

        //HYPRE_StructMatrixSetStencil(A, stencil);
        HYPRE_StructMatrixInitialize(A);
    }*/

    /*Return an updated flux limited value for the diffusion coefficient */
    //use the gradient of the energy
//the energy
//the total opacity
// int ifreq frequency bin number
    double RadSolve::larsendelimiter(const Mesh &mesh, State &state, Pars &pars, double opact, int i, int j, int k, int ifreq,int ord)
    {
        double dif=0;
    
        //double D = 1.0 / (3.0 * opact); // Corrected diffusion coefficient
        double grad_magnitude, egradt=0;
        //the grad here is a pointer to the gradient vector
       
        int ic = index(i, j, k,pars);
        grad_magnitude = std::sqrt(
            grad_energy[ifreq][ic][0] * grad_energy[ifreq][ic][0] +
            grad_energy[ifreq][ic][1] * grad_energy[ifreq][ic][1] +
            grad_energy[ifreq][ic][2] * grad_energy[ifreq][ic][2]
        );
    
        egradt=grad_magnitude/(state.radiation_flux[ifreq][ic]+1e-20); // avoid div by 0

        if(ord ==2)
            return(std::sqrt(1.0/(((9.0*opact*opact) + egradt*egradt))));
        else
            return(1.0/std::pow((std::pow(3.0*opact,ord) + std::pow(egradt,ord)),1/ord));

        return dif;
    }

    //Gradient of the energy
    double RadSolve::divergence(const Mesh &mesh, State &state, Pars &pars, int i, int j, int k, int ifreq)
    {
        double div=0;
        int ic1,ic2;

        ic2 = index(i, j, k,pars);


        //compute the divergence of the (energy gradient multiplied by the corrected diffusion coefficient)
        for(int di=-1; di<=1; di+=2) {
            if(i+di>=0 && i+di<pars.nx)
                {
                    ic1 = index(i+di, j, k,pars);                    
                    div+= (diff_coeff[ifreq][ic1]*grad_energy[ifreq][ic1][0] - diff_coeff[ifreq][ic2]*grad_energy[ifreq][ic2][0])/(2.0*pars.dx);
                }        
        }

        
        for(int dj=-1; dj<=1; dj+=2) {
            if(j+dj>=0 && j+dj<pars.ny)
            {
                ic1 = index(i, j+dj, k,pars);
                div+= (diff_coeff[ifreq][ic1]*grad_energy[ifreq][ic1][1] - diff_coeff[ifreq][ic2]*grad_energy[ifreq][ic2][1])/(2.0*pars.dy);
            }
        }

        for(int dk=-1; dk<=1; dk+=2) {
            if(k+dk>=0 && k+dk<pars.nz)
            {
                ic1 = index(i, j, k+dk,pars);
                div+= (diff_coeff[ifreq][ic1]*grad_energy[ifreq][ic1][2] - diff_coeff[ifreq][ic2]*grad_energy[ifreq][ic2][2])/(2.0*pars.dz); 
            }       
        }
        
        return div; // A small value to avoid division by zero

    }

    double RadSolve::gradenergy(const Mesh& mesh, State& state, Pars &pars)
    {
        //compute gradenergy for each group in each cell
        // Gradient of energy field  // [group][cell][3] 
        //std::vector<std::vector<std::vector<double>>> grad_energy

        double avg=0; // A small value to avoid division by zero
        double grad1=0,grad2=0,grad3=0;
    
        
        int ic;
        for(int n=0; n<pars.num_freq_bins; n++) {
            for (int k = 0; k < pars.nz; ++k)
                for (int j = 0; j < pars.ny; ++j)
                    for (int i = 0; i < pars.nx; ++i) {

                        ic=index((i+1<pars.nx?i+1:i),j,k,pars);
                        grad1=(state.radiation_flux[n][ic]-state.radiation_flux[n][ic])/(2*pars.dx);
                        ic=index(i,(j+1<pars.ny?j+1:j),k,pars);
                        grad2=(state.radiation_flux[n][ic]-state.radiation_flux[n][ic])/(2*pars.dy);
                        ic=index(i,j,(k+1<pars.nz?k+1:k),pars);
                        grad3=(state.radiation_flux[n][ic]-state.radiation_flux[n][ic])/(2*pars.dz);
                        avg+=grad1+grad2+grad3;
                        grad_energy[n][ic][0]=grad1;
                        grad_energy[n][ic][1]=grad2;
                        grad_energy[n][ic][2]=grad3;
        }
    }

    return avg;

    }


    void RadSolve::solveRadiationTransport(const Mesh& mesh, State& state, Pars &pars, double t) {


        std::mt19937 g(rd());
        std::vector<int> shuffled = ordered;
        std::shuffle(shuffled.begin(), shuffled.end(), g);
       
        std::cout<<"In solveRadiationTransport"<<std::endl;
        gradenergy(mesh,state,pars); // Compute the gradient of the energy
        std::cout<<"Computed gradenergy"<<std::endl;
        std::vector<double> values(7);
        std::vector<double> E_new(mesh.num_cells, 0.0);
       
        double a_nu = 0.7; // Initialize a_nu  should be computed using rosseland mean opacity
        double kappa_nu ;
        double sumrhs=0.0;
        double sumcdt=0.0;
        double sumd=0.0;
        double summu=0.0;
        double sumdiag=0.0;
        double sumemis=0.0;
        int idx=0;
        double sum1=0.0,sum2=0.0,sum3=0.0,sum4=0.0;

        int iindex[7]={0,1,2,3,4,5,6};
        int ilower[3] = {0, 0, 0};
        int iupper[3] = {pars.nx - 1, pars.ny - 1, pars.nz - 1};
        int offsets[7][3] = {{0,0,0},{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};


        for(int nf=0; nf<pars.num_freq_bins; nf++) {
            int n= shuffled[nf] - 1; // Get the shuffled frequency index (0-based)

            std::cout<<"Solving frequency bin "<<n<<std::endl;
        for (int k = 0; k < pars.nz; ++k)
        for (int j = 0; j < pars.ny; ++j)
        for (int i = 0; i < pars.nx; ++i) {
                                        idx = index(i, j, k,pars);
                                        sum1=0;
                                        sum2=state.heat_capacity[idx];
                                        sum3=0.0;
                                        sum4=0.0;
                                        ddelr[n][idx]=0.0;
                                        diff_coeff[n][idx]=0.0;
                                        for(int nf1=0; nf1<pars.num_freq_bins; nf1++) {
                                            sum1+=c*state.sigma_a[nf1][idx]*state.dBnudT(state.temperature[idx],(nf1+1)*1.0e14);
                                            sum3=c*state.sigma_a[nf1][idx]*(state.Bag[nf1][idx]);
                                            sum4=sum4+c*state.sigma_a[nf1][idx]*state.radiation_flux[nf1][idx];

                                        }
                                        diff_coeff[n][idx]=larsendelimiter(mesh,state,pars,state.sigma_a[n][idx],i,j,k,n);
                                        ddelr[n][idx]=divergence(mesh,state,pars,i,j,k,n);
                                        ;//std::cout<<"Computed diff_coeff and ddelr "<<diff_coeff[n][idx]<<" "<<ddelr[n][idx]<<std::endl;
                                        kappa_nu=1.0/(sum2+sum1);

                                        double diag = pars.scale*((1.0/ (c*pars.dt)) + ddelr[n][idx] + state.sigma_a[n][idx]);
                                        double D=ddelr[n][idx];
                                        values[0] = diag;
                                        for (int s = 1; s < 7; ++s) values[s] = -D ;

                                        int ijk[3] = {i, j, k};
                                        //std::cout<<"Setting matrix and emission values at cell "<<idx<<" i,j,k "<<i<<" "<<j<<" "<<k<<std::endl
                                         HYPRE_StructMatrixSetValues(A, ijk, 7, (HYPRE_Int*)iindex, values.data());
                                        double emission = pars.emisscale*state.sigma_a[n][idx] * state.Bag[n][idx]; // if treating B_nu as external source
                                        ;//std::cout<<"Computed emission "<<emission<<std::endl;
                                        //emission=0.0;
                                        double rhs = (state.radiation_flux[n][idx]  / (c*pars.dt)) + state.sigma_a[n][idx]*state.Bag[n][idx] +emission;

                                        //std::cout<< rhs<<" "<<emission<<" "<<diag<<" "<<D<<" "<<state.sigma_a[n][idx]<<" "<<state.radiation_flux[n][idx]<<" "<<Bag[n][idx]<<" "<<pars.dt<<std::endl;

                                        sumrhs=(fabs(rhs)>sumrhs?fabs(rhs):sumrhs);
                                        sumemis=(fabs(emission)>sumemis?fabs(emission):sumemis);
                                        sumdiag=(fabs(diag)>sumdiag?fabs(diag):sumdiag);
                                        sumcdt=(fabs(pars.scale*(1.0/ (c*pars.dt)))>sumcdt?fabs(pars.scale*(1.0/ (c*pars.dt))):sumcdt);
                                        summu=(fabs(pars.scale*state.sigma_a[n][idx])>summu?fabs(pars.scale*state.sigma_a[n][idx]):summu);
                                        sumd=(fabs(pars.scale*D/(pars.dx*pars.dx))>sumd?fabs(pars.scale*D/(pars.dx*pars.dx)):sumd);
                                        //HYPRE_StructVectorSetValues(b, ijk, &rhs);
                                        HYPRE_StructVectorSetValues(b, ijk, rhs);

                    }  //loop over cells

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

                    for (int k = 0; k < pars.nz; ++k)
                        for (int j = 0; j < pars.ny; ++j)
                            for (int i = 0; i < pars.nx; ++i) {
                                idx = index(i, j, k,pars);
                                state.radiation_fluxn[n][idx] = E_new[idx];
                            }   


        }//outer loop over frequencies   

    }

    void RadSolve::apply_milne_boundary_conditions(Mesh& mesh, State& state, Pars &pars)
    {

        int i,j,k;
        double D;
        double emis;
        double F_in, F_out, refg;
        refg=pars.refg; //reflection factor 0=absorbing, 1=perfect reflection

        //here we'll compute lates D and emission values for the boundary cells
        // Apply Milne boundary conditions at the domain boundaries
        for (const auto& bc : mesh.boundaries) {
            int cell = bc.cell;
            if (bc.type == INLET) {
                // For inlet, set radiation flux to a fixed value (e.g., blackbody at TINI)
                for (int g = 0; g < NUM_GROUPS; ++g) {
                    state.radiation_flux[g][cell] = STEFAN_BOLTZMANN * std::pow(TINI, 4) / NUM_GROUPS; // Simplified
                }
            } else if (bc.type == OUTLET) {
                // For outlet, apply zero-gradient condition
                for (int g = 0; g < NUM_GROUPS; ++g) {
                    // Assuming a simple 1D layout for illustration
                    int neighbor = cell + 1; // This should be determined based on mesh connectivity
                    if (neighbor < mesh.num_cells) {
                        state.radiation_flux[g][cell] = state.radiation_flux[g][neighbor];
                    }
                }
            } else if (bc.type == WALL) {
                // For wall, set radiation flux to zero or reflective condition
                i=cell%pars.nx;
                j=cell/pars.nx;
                k=cell/(pars.nx*pars.ny);  
                //be careful and check the signs here for emission and incoming flux
                for (int g = 0; g < pars.num_freq_bins; ++g) {              
                    emis=pars.emisscale*B_nu(state.temperature[cell],(g+1)*pars.df); // if treating B_nu as external source
                    D=larsendelimiter(mesh,state,pars,state.sigma_a[g][cell],i,j,k,g);
                    //state.radiation_flux[g][cell] = 0.0; // Absrbing wall
                    if(bc.wall_type==WUY) // wall
                    {
                        F_in=D*(grad_energy[g][cell][1]); // incoming flux note these are the old energy gradients
                        F_out=emis + (refg)*F_in; // outgoing flux
                        state.radiation_flux[g][cell] += F_out*pars.dy/D;
                    }                        
                    else if(bc.wall_type==WLY) //lower wall
                    {
                        F_in=-D*(grad_energy[g][cell][1]); // incoming flux note these are the old energy gradients
                        F_out=emis + (refg)*F_in; // outgoing flux
                        state.radiation_flux[g][cell] += F_out*pars.dy/D;
                    }                    
                    else if(bc.wall_type==WUX) //lower wall
                    {
                        F_in=D*(grad_energy[g][cell][0]); // incoming flux note these are the old energy gradients
                        F_out=emis + (refg)*F_in; // outgoing flux
                        state.radiation_flux[g][cell] += F_out*pars.dx/D;
                    }                        
                    else if(bc.wall_type==WLX) //upper wall
                    {
                        F_in=-D*(grad_energy[g][cell][0]); // incoming flux note these are the old energy gradients
                        F_out=emis + (refg)*F_in; // outgoing flux
                        state.radiation_flux[g][cell] += F_out*pars.dx/D;
                    }
                    else
                        state.radiation_flux[g][cell] = emis; // default to emissive only
                }
            }
        }
    }


    void RadSolve::apply_reflect_boundary_conditions(Mesh& mesh, State& state, Pars &pars)
    {

        // Apply Milne boundary conditions at the domain boundaries
        for (const auto& bc : mesh.boundaries) {
            int cell = bc.cell;
            if (bc.type == INLET) {
                // For inlet, set radiation flux to a fixed value (e.g., blackbody at TINI)
                for (int g = 0; g < NUM_GROUPS; ++g) {
                    state.radiation_flux[g][cell] = STEFAN_BOLTZMANN * std::pow(TINI, 4) / NUM_GROUPS; // Simplified
                }
            } else if (bc.type == OUTLET) {
                // For outlet, apply zero-gradient condition
                for (int g = 0; g < NUM_GROUPS; ++g) {
                    // Assuming a simple 1D layout for illustration
                    int neighbor = cell + 1; // This should be determined based on mesh connectivity
                    if (neighbor < mesh.num_cells) {
                        state.radiation_flux[g][cell] = state.radiation_flux[g][neighbor];
                    }
                }
            } else if (bc.type == WALL) {
                // For wall, set radiation flux to zero or reflective condition
                for (int g = 0; g < NUM_GROUPS; ++g) {
                    state.radiation_flux[g][cell] = 0.0; // Absorbing wall
                }
            }
        }
    }
