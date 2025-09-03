

#pragma once

constexpr int NX = 160, NY = 64, NZ = 1;
constexpr int NSTEP = 50000;
constexpr int N_SAVEINTERVAL=50;
const double DX = 5*0.07 / NX; //these units are in cm  should be 7cm
const double DY = 5*0.05 / NY;  //should be 5cm
const double DZ = 1.0 / NZ;  //should be 0.5cm
constexpr double dt = 1.0e-11;//1.0e-11;  // time step size
const int num_freq_bins = 10; // Number of frequency bins
const int NUM_GROUPS = num_freq_bins; // Number of frequency bins
const double c = 3e10; // Speed of light in cm/s
const double a = 7.5646e-15; // Radiation constant in erg/cm^3/K^4
const double sigma = 5.6704e-5; // Stefan-Boltzmann constant in erg/cm^2/s/K^4
const double STEFAN_BOLTZMANN=sigma;
const double h = 6.626e-27; // Planck's constant in erg.s
const double k = 1.38064852e-16; // Boltzmann constant in erg/K
const double TMIN = 293.0; // Minimum temperature - comfortable room temperature
const double TMAX = 10000.0; // Maximum temperature
const double TINI = 100000.0;//    //1.16e7; // Kelvin   11604525.0061598; // Maximum temperature
const double EINILT2 = 0.000000000001; // Initial temperature for top hot configuration
const double EINIEQ2 = 0.0000000001; // Initial temperature for equilibrium configuration
const int    BOUNDTYPE = 1; // Boundary type for the fixed temp is 1 reflected energy i1 2, and background is 0, 3 do nothing
const double temptol=0.1;
const int maxtiter=10;  //temperature iteration
const int BUY=1;//upper y 1
const int BLY=2;//lower y 2
const int BLX=3;//left x 3
const int BRX=4;//right x 4
const double SCALE=1.0;// Scaling factor for the diagonal term in the matrix
const double EMISSCALE=1.0; // Scaling factor for the emission term


//solver
#define GNX 50  // Grid points in X
#define GNY 50  // Grid points in Y
#define TIME_STEPS 10 // Number of time steps
#define DT 0.01 // Time step size
#define ALPHA 0.1 // Absorption coefficient

//physics state
constexpr double STEFAN_BOLTZMANN = 5.670374419e-8; // W·m⁻²·K⁻⁴
constexpr int NUM_GROUPS = 3; // Example: 3 frequency groups
constexpr double DT = 1e-9;   // Time step in seconds