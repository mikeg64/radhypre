#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <HYPRE.h>
#include <HYPRE_struct_ls.h>
#include "ex.h"
#include "setup.h"
#include "geometry.h"   //defines mesh
#include "material.h"   //defines material properties
#include "physics.h"



void solve_radiation_groups(const Mesh& mesh, State& state);
int setupsolver() ;

class RadSolve {
public:
    RadSolve(int mnx, int mny);
    ~RadSolve();
    //const std::vector<std::vector<double>>& getTemperatureField() const;
    //void setTemperature(int i, int j, double value);
    //void readMesh(const std::string& filename);
    //void setupGrid();
    void solveRadiationTransport(const Mesh& mesh, State& state, double t);






private:
    int nx, ny; // Grid dimensions 
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    std::vector<std::vector<double>> T; // Temperature field
    HYPRE_StructSolver solver;
    std::vector<std::pair<double, double>> nodes;
};
