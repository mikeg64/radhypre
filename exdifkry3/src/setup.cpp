#include "../include/setup.h"


// Simulation parameters
Pars::Pars() {
   
    // Default parameters
    nx = NX;
    ny = NY;
    nz = NZ;
    nstep = NSTEP;
    nsaveinterval = N_SAVEINTERVAL;
    dx = DX;
    dy = DY;
    dz = DZ;
    dt = DT;
    df = 1.0e14; // frequency bin size
    refg=0.6; //reflection factor for walls 0=absorbing, 1=perfect reflection
    dtmin=DTMIN;
    dtmax=DTMAX; 
    num_freq_bins = NUM_FREQ_BINS;
    tmin = TMIN;
    tmax = TMAX;
    tini = TINI;
    einilt2 = EINILT2;
    einieq2 = EINIEQ2;
    boundtype = BOUNDTYPE;
    temptol = TEMPTOL;
    maxtiter = MAXITER;
    buy = BUY;
    bly = BLY;
    blx = BLX;
    brx = BRX;
    scale = SCALE;
    emisscale = EMISSCALE;
}

void Pars::readpars(std::string filename){
    std::ifstream infile(filename);
    if (!infile) {
        std::cout << "Error opening parameter file: " << filename << std::endl;
        return;
    }
    std::string line;
    while (std::getline(infile, line)) {    
        std::istringstream iss(line);
        std::string key;
        if (iss >> key) {
            if (key == "NX") iss >> nx;
            else if (key == "NY") iss >> ny;
            else if (key == "NSTEP") iss >> nstep;
            else if (key == "N_SAVEINTERVAL") iss >> nsaveinterval;
            else if (key == "DX") iss >> dx;
            else if (key == "DY") iss >> dy;
            else if (key == "DZ") iss >> dz;
            else if (key == "dt") iss >> dt;
            else if (key == "dtmin") iss >> dtmin;
            else if (key == "dtmax") iss >> dtmax;  
            else if (key == "num_freq_bins") iss >> num_freq_bins;
            else if (key == "TMIN") iss >> tmin;
            else if (key == "TMAX") iss >> tmax;
            else if (key == "TINI") iss >> tini;
            else if (key == "EINILT2") iss >> einilt2;
            else if (key == "EINIEQ2") iss >> einieq2;
            else if (key == "BOUNDTYPE") iss >> boundtype;
            else if (key == "temptol") iss >> temptol;
            else if (key == "maxtiter") iss >> maxtiter;
            else if (key == "BUY") iss >> buy;
            else if (key == "BLY") iss >> bly;
            else if (key == "BLX") iss >> blx;
            else if (key == "BRX") iss >> brx;
            else if (key == "SCALE") iss >> scale;
            else if (key == "EMISSCALE") iss >> emisscale;
        }
    }
}
