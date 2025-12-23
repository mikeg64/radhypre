#pragma once

#include <vector>
#include "setup.h"





class Cell {
    public:
        double x, y;
        bool in_pipe;  // true if cell is inside the pipe to be replaced by region also has integere to denote which region it's in see mesh
        int material_id;
};

 

enum BoundaryType {
    WALL,
    INLET,
    OUTLET
};

enum WallType {
    WUY,
    WLY,
    WUX,
    WLX
};


 

class BoundaryCondition {
public:
    int cell;
    BoundaryType type;
    WallType wall_type;
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
    std::vector<Cell> cells;  /// list of cells in the mesh
                              //also have identifier for list of regions
    std::vector<BoundaryCondition> boundaries;

};

Mesh setup_crooked_pipe_geometry(Pars &pars);
bool is_wall_cell(const Mesh& mesh, int idx); 




