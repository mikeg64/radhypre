#pragma once

#include "setup.h"

#include "geometry.h"   //defines mesh

// Physical constants

constexpr double STEFAN_BOLTZMANN = 5.670374419e-8; // W·m⁻²·K⁻⁴

constexpr int NUM_GROUPS = 3; // Example: 3 frequency groups

constexpr double DT = 1e-9;   // Time step in seconds

struct PhysicsState {
    std::vector<std::vector> radiation_flux;   // [group][cell]</std::vector
    std::vector<std::vector> sigma_a;          // [group][cell]</std::vector
    std::vector<std::vector> source_term;      // [group][cell]</std::vector
    std::vector temperature;                   // [cell]
    std::vector heat_capacity;                 // [cell]
};


void linearize_emissive_source(PhysicsState& state)
void solve_material_heating(const Mesh& mesh, PhysicsState& state);

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
double B_nu(double T, double nu, double dbnu=-1.0) {
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
double larsendelimiter(double opact, double Eg[num_freq_bins][NX][NY][NZ],double GradEg[num_freq_bins][NX][NY][NZ][3], int i, int j, int k, int ifreq,int ord=2) {
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
double gradenergy(double GradEg[num_freq_bins][NX][NY][NZ][3], double Eg[num_freq_bins][NX][NY][NZ]) {
    double avg=0; // A small value to avoid division by zero
    double grad1=0,grad2=0,grad3=0;
    

    for(int n=0; n<num_freq_bins; n++) {
        for (int k = 0; k < NZ; ++k)
        for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i) {
            grad1=(Eg[n][(i+1<NX?i+1:i)][j][k]-Eg[n][(i>0?i-1:i)][j][k])/(2*dx);
            grad2=(Eg[n][i][(j+1<NY?j+1:j)][k]-Eg[n][i][(j>0?j-1:j)][k])/(2*dy);
            grad3=(Eg[n][i][j][(k+1<NZ?k+1:k)]-Eg[n][i][j][(k>0?k-1:k)])/(2*dz);
            avg+=grad1+grad2+grad3;
            GradEg[n][i][j][k][0]=grad1;
            GradEg[n][i][j][k][1]=grad2;
            GradEg[n][i][j][k][2]=grad3;
        }
    }
    return avg/(3*NX*NY*NZ*num_freq_bins); // A small value to avoid division by zero
}


//Gradient of the energy
double divergence(double diffcoeff[num_freq_bins][NX][NY][NZ], double gradeg[num_freq_bins][NX][NY][NZ][3],int i, int j, int k, int ifreq) {
   
    double grad=0;
    //compute the divergence of the (energy gradient multiplied by the corrected diffusion coefficient)
    for(int di=-1; di<=1; di+=2) {
        if(i+di>=0 && i+di<NX)
            grad+= (diffcoeff[ifreq][i+di][j][k]*gradeg[ifreq][i+di][j][k][0] - diffcoeff[ifreq][i][j][k]*gradeg[ifreq][i][j][k][0])/(2.0*dx);        
    }

    
    for(int dj=-1; dj<=1; dj+=2) {
        if(j+dj>=0 && j+dj<NY)
            grad+= (diffcoeff[ifreq][i][j+dj][k]*gradeg[ifreq][i][j+dj][k][1] - diffcoeff[ifreq][i][j][k]*gradeg[ifreq][i][j][k][1])/(2.0*dy);
    }

    for(int dk=-1; dk<=1; dk+=2) {
        if(k+dk>=0 && k+dk<NZ)
            grad+= (diffcoeff[ifreq][i][j][k+dk]*gradeg[ifreq][i][j][k+dk][2] - diffcoeff[ifreq][i][j][k]*gradeg[ifreq][i][j][k][2])/(2.0*dz);        
    }
    
    return grad; // A small value to avoid division by zero
}