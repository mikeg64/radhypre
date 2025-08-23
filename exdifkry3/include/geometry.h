#pragma once

#include <vector>
#include "setup.h"

 

struct Cell {

    double x, y;

    bool in_pipe;

    int material_id;

};

 

enum BoundaryType {

    WALL,

    INLET,

    OUTLET

};

 

struct BoundaryCondition {

    int cell;

    BoundaryType type;

};

 

struct GMesh {

    int nx, ny;

    double dx, dy;

    int num_cells;

    std::vector<Cell> cells;

    std::vector<BoundaryCondition> boundaries;

};


class Mesh {


public:
    Mesh(int nx, int ny);
    ~Mesh();
    const std::vector<std::vector<double>>& getTemperatureField() const;
    void setTemperature(int i, int j, double value);






private:

    int nx, ny;
    double dx, dy;
    int num_cells;
    std::vector<Cell> cells;
    std::vector<BoundaryCondition> boundaries;

};