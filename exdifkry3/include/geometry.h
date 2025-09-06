#pragma once

#include <vector>
#include "setup.h"





class Cell {
    public:
        double x, y;
        bool in_pipe;
        int material_id;
};

 

enum BoundaryType {
    WALL,
    INLET,
    OUTLET
};

 

class BoundaryCondition {
public:
    int cell;
    BoundaryType type;
};

class Mesh {


public:
    Mesh(int mnx, int mny, double dx, double dy);
    ~Mesh();
    //const std::vector<std::vector<double>>& getTemperatureField() const;
   // void setTemperature(int i, int j, double value);






public:

    int nx, ny;
    double dx, dy;
    int num_cells;
    std::vector<Cell> cells;
    std::vector<BoundaryCondition> boundaries;

};

Mesh setup_crooked_pipe_geometry();
bool is_wall_cell(const Mesh& mesh, int idx); 




