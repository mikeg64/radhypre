#include "../include/physics.h"

 

void linearize_emissive_source(const Mesh& mesh, State& state, Pars &pars) {
    for (int i = 0; i < mesh.num_cells; ++i) {
        double T = state.temperature[i];
        double T4 = std::pow(T, 4);
        for (int g = 0; g < NUM_GROUPS; ++g) {
            state.source_term[g][i] = state.sigma_a[g][i] * STEFAN_BOLTZMANN * T4;
        }
    }
}

 

void solve_material_heating(const Mesh& mesh, State& state, Pars &pars) {


        double etotn=0.0;
        double sum1, sum2, sum4;
        
        //compute new temps
        double delTmax=0.0;
        double delT;
        for (int k = 0; k < pars.nz; ++k)
        for (int j = 0; j < pars.ny; ++j)
        for (int i = 0; i < pars.nx; ++i) {        

                sum1=0;
                sum4=state.heat_capacity[index(i,j,k,pars)];
                sum2=0.0;
                int ic=index(i,j,k,pars);
                for(int n=0; n<pars.num_freq_bins; n++) 
                {
                    sum1+= (c*state.sigma_a[n][ic]*(state.radiation_flux[n][ic]-state.Bag[n][ic]));
                    sum2+= (c*state.sigma_a[n][ic]*dBnudT(state.temperature[ic],(n+1)*1.0e14));
                    etotn+=state.radiation_flux[n][index(i,j,k,pars)];
                              
                }
                //Tn[i][j][k]=Tc[i][j][k] + ((sum1/sum4)/(1+sum2));
                state.temperature[ic] += pars.dt*(sum1/(sum2+sum4));
                delT=fabs(pars.dt*(sum1/(sum2+sum4)));
                if(delT>delTmax) delTmax=delT;
                //apply boundary condition
                etotn += state.heat_capacity[index(i,j,k)]*state.temperature[ic]; // Calculate total energy at current time step
            }
            state.etot=etotn;
            std::cout<<"max temp change "<<delTmax<<std::endl;



}

/*void solve_radiation_groups(const Mesh& mesh, State& state) {
    // Placeholder for radiation transport solver
    // This would typically involve solving a system of equations for each group
    for (int g = 0; g < NUM_GROUPS; ++g) {
        for (int i = 0; i < mesh.num_cells; ++i) {
            // Simple diffusion approximation (placeholder)
            state.radiation_flux[g][i] += DT * (-state.sigma_a[g][i] * state.radiation_flux[g][i] + state.source_term[g][i]);
        }
    }
}*/

State initialize_physics(Mesh& mesh,  Materials& materials) {

    State state(mesh.nx, mesh.ny);
    int num_cells = mesh.num_cells;

    // Resize group-dependent vectors
    state.radiation_flux.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
    state.sigma_a.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
    state.source_term.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
    state.Bag.resize(NUM_GROUPS, std::vector(num_cells, 0.0));

    // Resize temperature and heat capacity
    state.temperature.resize(num_cells, 0.0);
    state.heat_capacity.resize(num_cells, 0.0);

    for (int i = 0; i < num_cells; ++i) {
        int mat_id = mesh.cells[i].material_id;
        // Set initial temperature (e.g., 1 keV blackbody)
        state.temperature[i] = TMIN; //1.16e7; // Kelvin
        // Set heat capacity from material database
        state.heat_capacity[i] = materials.get_heat_capacity(mat_id);
        // Set absorption coefficients per group
        for (int g = 0; g < NUM_GROUPS; ++g) {
            state.sigma_a[g][i] = materials.get_sigma_a(mat_id);// this for different frequency groups (mat_id, g);
        }

        // Initialize source term using blackbody emission
        //double T4 = std::pow(state.temperature[i], 4);
        state.etot=0.0;
        for (int g = 0; g < NUM_GROUPS; ++g) {
            //state.source_term[g][i] = state.sigma_a[g][i] * STEFAN_BOLTZMANN * T4;
            state.source_term[g][i] = planck_emission(  (g+1)*1e14,state.temperature[i]) ; 
            state.Bag[g][i] = planck_emission(  (g+1)*1e14,state.temperature[i]) ; 
            state.radiation_flux[g][i] = planck_emission(  (g+1)*1e14,state.temperature[i]) ; 
            state.etot+=state.radiation_flux[g][i];
        }
    }

    return state;

}


void apply_milne_boundary_conditions(Mesh& mesh, State& state)
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


void apply_reflect_boundary_conditions(Mesh& mesh, State& state)
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

    State::State(int nx, int ny)
    {
        int num_cells = nx * ny;
        // Resize group-dependent vectors
        radiation_flux.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
        radiation_fluxn.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
        sigma_a.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
        source_term.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
        Bag.resize(NUM_GROUPS, std::vector(num_cells, 0.0));

        // Resize temperature and heat capacity
        temperature.resize(num_cells, 0.0);
        heat_capacity.resize(num_cells, 0.0);

        etot=0.0;
    }

    State::State(State &other) : radiation_flux(other.radiation_flux), sigma_a(other.sigma_a), source_term(other.source_term), temperature(other.temperature), heat_capacity(other.heat_capacity) 
    {

    }


    State::~State(){}
    std::vector<double>& State::getTemperatureField()  {
        
        return temperature;
    }

    void State::setTemperature(int i, int j, double value){}


    void State::copy(const State& other)
    {
        temperature = other.temperature;
        heat_capacity = other.heat_capacity;
        radiation_flux = other.radiation_flux;
        radiation_fluxn = other.radiation_fluxn;
        sigma_a = other.sigma_a;
        source_term = other.source_term;
    }

    int State::getNumCells() const
    {
        return temperature.size();
    }

    int State::getNumGroups() const
    {
        return radiation_flux.size();
    }

    double State::getRadiationFlux(int group, int cell) const
    {
        return radiation_flux[group][cell];
    }


    void State::setRadiationFlux(int group, int cell, double value)
    {
        radiation_flux[group][cell] = value;
    }

    double State::getSigmaA(int group, int cell) const
    {
        return sigma_a[group][cell];
    }

    void State::setSigmaA(int group, int cell, double value)
    {
        sigma_a[group][cell] = value;
    }

    double State::getSourceTerm(int group, int cell) const
    {
        return source_term[group][cell];
    }

    void State::setSourceTerm(int group, int cell, double value)
    {
        source_term[group][cell] = value;
    }

    double State::getHeatCapacity(int cell) const
    {
        return heat_capacity[cell];
    }


    void State::setHeatCapacity(int cell, double value)
    {
        heat_capacity[cell] = value;
    }

   // Function to write temperature data to a VTK file
    void State::write_vtk_file(const std::vector<double> &scalf, int time_step, Pars &pars, const std::string &scalfieldname,const std::string &filename)
     {
        std::string vtk_filename = filename + "_" + std::to_string(time_step) + ".vtk";
        std::ofstream vtk_file(vtk_filename);

        if (!vtk_file) {
            std::cerr << "Error: Unable to open file for writing!" << std::endl;
            return;
        }

        vtk_file << "# vtk DataFile Version 3.0\n";
        vtk_file << scalfieldname << " field output\n";
        vtk_file << "ASCII\n";
        vtk_file << "DATASET STRUCTURED_POINTS\n";
        vtk_file << "DIMENSIONS " << pars.nx << " " << pars.ny << " 1\n";
        vtk_file << "ORIGIN 0 0 0\n";
        vtk_file << "SPACING 1.0 1.0 1.0\n";
        vtk_file << "POINT_DATA " << pars.nx * pars.ny << "\n";
        vtk_file << "SCALARS "<< scalfieldname <<" float\n";
        vtk_file << "LOOKUP_TABLE default\n";

        // Write the temperature values in row-major order
        int k=0;
        for (int j = 0; j < pars.ny; j++) {
            for (int i = 0; i < pars.nx; i++) {
                vtk_file << scalf[index(i,j,k,pars)] << "\n";
            }
        }

        vtk_file.close();
        std::cout << "VTK file '" << vtk_filename << "' written successfully!" << std::endl;
    }



    // Function to compute Planck distribution
    double State::dBnudT(double T, double nu)
    {
        return (exp(h*nu/(k*T))*2 * pow(h,2) * pow(nu, 4) / (k*pow(T,2)*pow(c, 2))) / pow(exp(h * nu / (k * T)) - 1,2);
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
double B_nu(double T, double nu, double dbnu) {
    double tot=0;
    double nud=nu;

    if(dbnu==-1.0) {
        tot = (2 * h * std::pow(nu, 3) / std::pow(c, 2)) / (std::exp(h * nu / (k * T)) - 1);
    } else {
        nud=nu-(dbnu/2.0);
        nu=(nud>0?nud:nu); // Avoid division by zero
        tot = (2 * h * std::pow(nu, 3) / std::pow(c, 2)) / (std::exp(h * nu / (k * T)) - 1) * dbnu;
    }


    tot= (2 * h * pow(nu, 3) / pow(c, 2)) / (exp(h * nu / (k * T)) - 1);

    return tot;
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


//return the corrected diffusion coefficient
//use the gradient of the energy
//the energy
//the total opacity
// int ifreq frequency bin number
double larsendelimiter(double opact, double Eg[NUM_FREQ_BINS][NX][NY][NZ],double GradEg[NUM_FREQ_BINS][NX][NY][NZ][3], int i, int j, int k, int ifreq,int ord) {
    //double D = 1.0 / (3.0 * opact); // Corrected diffusion coefficient
    double grad_magnitude, egradt=0;
    //the grad here is a pointer to the gradient vector
    
    grad_magnitude = std::sqrt(
        GradEg[ifreq][i][j][k][0] * GradEg[ifreq][i][j][k][0] +
        GradEg[ifreq][i][j][k][1] * GradEg[ifreq][i][j][k][1] +
        GradEg[ifreq][i][j][k][2] * GradEg[ifreq][i][j][k][2]
    );
    
    egradt=grad_magnitude/(Eg[ifreq][i][j][k]+1e-20); // avoid div by 0

    if(ord ==2)
        return(std::sqrt(1.0/(((9.0*opact*opact) + egradt*egradt))));
    else
        return(1.0/std::pow((std::pow(3.0*opact,ord) + std::pow(egradt,ord)),1/ord));


}


//Gradient of the energy
double gradenergy(double GradEg[NUM_FREQ_BINS][NX][NY][NZ][3], double Eg[NUM_FREQ_BINS][NX][NY][NZ]) {
    double avg=0; // A small value to avoid division by zero
    double grad1=0,grad2=0,grad3=0;
    

    for(int n=0; n<NUM_FREQ_BINS; n++) {
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {
            grad1=(Eg[n][(i+1<NX?i+1:i)][j][k]-Eg[n][(i>0?i-1:i)][j][k])/(2*DX);
            grad2=(Eg[n][i][(j+1<NY?j+1:j)][k]-Eg[n][i][(j>0?j-1:j)][k])/(2*DY);
            grad3=(Eg[n][i][j][(k+1<NZ?k+1:k)]-Eg[n][i][j][(k>0?k-1:k)])/(2*DZ);
            avg+=grad1+grad2+grad3;
            GradEg[n][i][j][k][0]=grad1;
            GradEg[n][i][j][k][1]=grad2;
            GradEg[n][i][j][k][2]=grad3;
        }
    }
    return avg/(3*NX*NY*NZ*NUM_FREQ_BINS); // A small value to avoid division by zero
}


//Gradient of the energy
double divergence(double diffcoeff[NUM_FREQ_BINS][NX][NY][NZ], double gradeg[NUM_FREQ_BINS][NX][NY][NZ][3],int i, int j, int k, int ifreq) {
   
    double grad=0;
    //compute the divergence of the (energy gradient multiplied by the corrected diffusion coefficient)
    for(int di=-1; di<=1; di+=2) {
        if(i+di>=0 && i+di<NX)
            grad+= (diffcoeff[ifreq][i+di][j][k]*gradeg[ifreq][i+di][j][k][0] - diffcoeff[ifreq][i][j][k]*gradeg[ifreq][i][j][k][0])/(2.0*DX);        
    }

    
    for(int dj=-1; dj<=1; dj+=2) {
        if(j+dj>=0 && j+dj<NY)
            grad+= (diffcoeff[ifreq][i][j+dj][k]*gradeg[ifreq][i][j+dj][k][1] - diffcoeff[ifreq][i][j][k]*gradeg[ifreq][i][j][k][1])/(2.0*DY);
    }

    for(int dk=-1; dk<=1; dk+=2) {
        if(k+dk>=0 && k+dk<NZ)
            grad+= (diffcoeff[ifreq][i][j][k+dk]*gradeg[ifreq][i][j][k+dk][2] - diffcoeff[ifreq][i][j][k]*gradeg[ifreq][i][j][k][2])/(2.0*DZ);        
    }
    
    return grad; // A small value to avoid division by zero
}

