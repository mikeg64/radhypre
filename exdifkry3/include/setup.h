

#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


// Simulation parameters
//nx is a factor of 7
//ny is a factor of 5
constexpr int NX = 140, NY = 60, NZ = 1;
const int NSTEP = 50000;
const int N_SAVEINTERVAL=50;
const double DX = 5*0.07 / NX; //these units are in cm  should be 7cm
const double DY = 5*0.05 / NY;  //should be 5cm
const double DZ = 1.0 / NZ;  //should be 0.5cm
const double DT = 1.0e-11;//1.0e-11;  // time step size
const double DTMIN=1.0e-12; //minimum time step
const double DTMAX=1.0e-7; //maximum time step  
constexpr int NUM_FREQ_BINS = 10; // Number of frequency bins
const int NUM_GROUPS = NUM_FREQ_BINS; // Number of frequency bins
const double c = 3e10; // Speed of light in cm/s
const double a = 7.5646e-15; // Radiation constant in erg/cm^3/K^4

//const double STEFAN_BOLTZMANN=sigma;
const double h = 6.626e-27; // Planck's constant in erg.s
const double k = 1.38064852e-16; // Boltzmann constant in erg/K
const double TMIN = 293.0; // Minimum temperature - comfortable room temperature
const double TMAX = 10000.0; // Maximum temperature
const double TINI = 10000000000.0;//    //1.16e7; // Kelvin   11604525.0061598; // Maximum temperature
const double EINILT2 = 0.000000000001; // Initial temperature for top hot configuration
const double EINIEQ2 = 0.0000000001; // Initial temperature for equilibrium configuration
const int    BOUNDTYPE = 1; // Boundary type for the fixed temp is 1 reflected energy i1 2, and background is 0, 3 do nothing
const double TEMPTOL=0.1;
const int MAXITER=10;  //temperature iteration
const int BUY=0;//upper y 1
const int BLY=1;//lower y 2
const int BLX=2;//left x 3
const int BRX=3;//right x 4
const double SCALE=1.0;// Scaling factor for the diagonal term in the matrix
const double EMISSCALE=1.0; // Scaling factor for the emission term


//solver
const int  GNX=50 ; // Grid points in X
const int  GNY=50;  // Grid points in Y
const int TIME_STEPS=10; // Number of time steps

const double  ALPHA=0.1; // Absorption coefficient

//physics state
const double STEFAN_BOLTZMANN = 5.670374419e-8; // W·m⁻²·K⁻⁴
//constexpr int NUM_GROUPS = 3; // Example: 3 frequency groups
//const double DT = dt;   // Time step in seconds
const double sigma = STEFAN_BOLTZMANN; // Stefan-Boltzmann constant in erg/cm^2/s/K^4

class Pars{

public:
    Pars();
    ~Pars() {}

    void readpars(std::string filename);
    

    int nx;
    int ny;
    int nz;
    double dx;
    double dy;
    double dz;
    double dt;
    double dtmin;
    double dtmax;
    double time;
    double df; // frequency bin size
    double refg; //reflection factor for walls 0=absorbing, 1=perfect reflection
    int nstep;
    int nsaveinterval;
    int num_freq_bins;
    int num_groups;
    double tmin; // Minimum temperature - comfortable room temperature
    double tmax; // Maximum temperature
    double tini;// Initial temperature
    double einilt2; // Initial temperature for top hot configuration
    double einieq2; // Initial temperature for equilibrium configuration
    int boundtype; // Boundary type for the fixed temp is 1 reflected energy i1 2, and background is 0, 3 do nothing
    double temptol;
    int maxtiter;  //temperature iteration
    int buy;//upper y 1
    int bly;//lower y 2
    int blx;//left x 3
    int brx;//right x 4
    double scale;// Scaling factor for the diagonal term in the matrix
    double emisscale; // Scaling factor for the emission term
};