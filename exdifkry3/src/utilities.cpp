

#include <cstdio>
#include "../include/utilities.h"
#include "../include/setup.h"





int index(int i, int j, int k) {
    return i + NX * (j + NY * k);
}

int index(int i, int j, int k, Pars &pars) {
    return i + pars.nx * (j + pars.ny * k);
}







void write_vtk_file(double T[NX][NY][NZ], int timestep, Pars &pars) {
    FILE *vtkFile;
    char filename[50];

    sprintf(filename, "radiation_output_%d.vtk", timestep);
    vtkFile = fopen(filename, "w");

    if (!vtkFile) {
        printf("Error: Could not open file for writing!\n");
        return;
    }

    // Write VTK header
    fprintf(vtkFile, "# vtk DataFile Version 3.0\n");
    fprintf(vtkFile, "Radiation Transport Simulation\n");
    fprintf(vtkFile, "ASCII\n");
    fprintf(vtkFile, "DATASET STRUCTURED_POINTS\n");
    fprintf(vtkFile, "DIMENSIONS %d %d %d\n", NX, NY, NZ);
    fprintf(vtkFile, "ORIGIN 0 0 0\n");
    fprintf(vtkFile, "SPACING 1 1 1\n");

    // Write data section
    fprintf(vtkFile, "POINT_DATA %d\n", (NX * NY*NZ));
    fprintf(vtkFile, "SCALARS Temperature float 1\n");
    fprintf(vtkFile, "LOOKUP_TABLE default\n");
 
    // Write temperature values in row-major order
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++)
        for (int k = 0; k < NZ; k++) {
            fprintf(vtkFile, "%f ", T[i][j][k]);
        }
        fprintf(vtkFile, "\n");
    }

    fclose(vtkFile);
    printf("VTK file written successfully: %s\n", filename);
}