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

 

class GMesh {
public:
    GMesh() : nx(40), ny(20), dx(1.0), dy(1.0) {
        num_cells = nx * ny;
        cells.resize(num_cells);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                int idx = j * nx + i;
                cells[idx].x = i * dx;
                cells[idx].y = j * dy;
                // Define crooked pipe: horizontal then vertical bend   ?????
                bool in_pipe = (j >= 8 && j <= 12 && i < 30) || (i >= 28 && i <= 32 && j >= 8 && j <= 18);
                cells[idx].in_pipe = in_pipe;
                cells[idx].material_id = in_pipe ? 1 : 0;
            }
        }
        // Define boundaries
        for (int i = 0; i < num_cells; ++i) {
            if (cells[i].in_pipe) {
                if (is_wall_cell(i)) {
                    BoundaryCondition bc;
                    bc.cell = i;
                    bc.type = WALL;
                    boundaries.push_back(bc);
                }
            }
        }
    }
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