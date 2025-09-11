#pragma once

#include "setup.h"

#include "geometry.h"   //defines mesh
#include "material.h"   //defines material properties
//#include "solver.h"     //defines solver class

#include <cmath>
#include <vector>
#include <functional>
#include <iostream>

// Physical constants




class State {


public:
    State(int nx, int ny);
    State(State &other) ;

    ~State();
    std::vector<double>& getTemperatureField();
    void setTemperature(int i, int j, double value);
    void copy(const State& other);
    int getNumCells() const;
    int getNumGroups() const;
    double getRadiationFlux(int group, int cell) const; 
    void setRadiationFlux(int group, int cell, double value);
    double getSigmaA(int group, int cell) const;    
    void setSigmaA(int group, int cell, double value);
    double getSourceTerm(int group, int cell) const;
    void setSourceTerm(int group, int cell, double value);
    double getHeatCapacity(int cell) const;
    void setHeatCapacity(int cell, double value);







public:   //set methods and make private
    

    std::vector<std::vector<double>> radiation_flux;   // [group][cell]
    std::vector<std::vector<double>> sigma_a;          // [group][cell]
    std::vector<std::vector<double>> source_term;      // [group][cell]
    std::vector<double> temperature;                   // [cell]
    std::vector<double> heat_capacity;                 // [cell]
};







struct PhysicsState {
    std::vector<std::vector<double>> radiation_flux;   // [group][cell]</std::vector
    std::vector<std::vector<double>> sigma_a;          // [group][cell]</std::vector
    std::vector<std::vector<double>> source_term;      // [group][cell]</std::vector
    std::vector<double> temperature;                   // [cell]
    std::vector<double> heat_capacity;                 // [cell]
};


void linearize_emissive_source(const Mesh& mesh,State& state);
void solve_material_heating(const Mesh& mesh, State& state);
//void solve_radiation_groups(const Mesh& mesh, State& state); // see solver class
State initialize_physics(Mesh& mesh,  Materials& materials);
//same as B_nu but for Planck's law
// This function computes the spectral radiance of a black body at frequency nu and temperature T


void apply_milne_boundary_conditions(Mesh& mesh, State& state);
void apply_reflect_boundary_conditions(Mesh& mesh, State& state);

double planck_emission(double nu, double T);


// Function to compute Planck distribution
double B_nu(double T, double nu, double dbnu=-1.0);

// Function to compute Planck distribution
double dBnudT(double T, double nu);


double dBdT(double nu, double T) ;


// Define the opacity function
std::function<double(double)> make_kappa_nu(double kappa0, double nu0, double alpha);
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
    int N  // number of integration points
) ;

double rosseland_mean_opacity_group(
    double kappa, // frequency-dependent opacity
    double T,
    int num_freq_bins  // number of integration points
) ;






double planck_source(double T) ;

//return the corrected diffusion coefficient
//use the gradient of the energy
//the energy
//the total opacity
// int ifreq frequency bin number
double larsendelimiter(double opact, double Eg[NUM_FREQ_BINS][NX][NY][NZ],double GradEg[NUM_FREQ_BINS][NX][NY][NZ][3], int i, int j, int k, int ifreq,int ord=2) ;

//Gradient of the energy
double gradenergy(double GradEg[NUM_FREQ_BINS][NX][NY][NZ][3], double Eg[NUM_FREQ_BINS][NX][NY][NZ]);


//Gradient of the energy
double divergence(double diffcoeff[NUM_FREQ_BINS][NX][NY][NZ], double gradeg[NUM_FREQ_BINS][NX][NY][NZ][3],int i, int j, int k, int ifreq) ;
