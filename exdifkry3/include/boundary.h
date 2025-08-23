#pragma once



#include "geometry.h"
#include "physics.h"



 

double compute_milne_reflected_flux(double incident_flux, double temperature, double sigma_a);
void apply_milne_boundary_conditions(const Mesh& mesh, PhysicsState& state);

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
