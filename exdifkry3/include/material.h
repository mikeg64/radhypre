#pragma once


// Setting up the material properties for the simulation


double absorption_coeff(int i, int j) {
    double siga;
    double sigaupper=2000.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        siga = sigaupper ;
    } else {
        siga = 20.0;  //20.0
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

double initial_temperature(int i, int j) {
    double tini=300;
    if(i<4 &&  j< 3*(NY/2)/4   && j>NY/4) {
    //if(i<4 &&  j< (NY/2)/2   ) {
        tini=TMIN+1.0e9*std::exp(-((i-2)*(i-2)+(j-(NY)/2)*(j-(NY)/2))/4)*TINI; // top hot configuration
        std::cout << " Initial hot temp " << tini << std::endl;
    } else {
        tini=300.0;
    }

    return tini;
}
