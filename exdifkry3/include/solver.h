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




int setupsolver() {
    // This function sets up the solver parameters and initializes necessary variables
    // It can be expanded to include more complex setup logic if needed
    // For now, it simply returns true to indicate successful setup
    return 0;
}


class RadSolve {
public:
    RadSolve(int mnx, int mny);
    ~RadSolve();
    const std::vector<std::vector<double>>& getTemperatureField() const;
    void setTemperature(int i, int j, double value);
    void readMesh(const std::string& filename);
    void setupGrid();
    void solveRadiationTransport();






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
