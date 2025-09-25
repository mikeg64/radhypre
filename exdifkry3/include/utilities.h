#pragma once

#include <cstdio>
#include <iostream>
#include <fstream>

#include "../include/setup.h"


int index(int i, int j, int k);
int index(int i, int j, int k, Pars &pars);


void write_vtk_file(double T[NX][NY][NZ], int timestep, Pars &pars);