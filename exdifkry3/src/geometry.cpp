#include "../include/geometry.h"





Mesh::~Mesh() {}

Mesh::Mesh(int mnx, int mny, double dx, double dy)  : nx(mnx), ny(mny), dx(dx), dy(dy)
{
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
            if (!cells[i].in_pipe) {
                
                    BoundaryCondition bc;
                    bc.cell = i;
                    bc.type = WALL;
                    boundaries.push_back(bc);
                
            }
        }
}




Mesh setup_crooked_pipe_geometry() {
    Mesh mesh(NX, NY, DX, DY);
    
    //mesh.num_cells = mesh.nx * mesh.ny;
    //mesh.cells.resize(mesh.num_cells);

 

    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            int idx = j * mesh.nx + i;
            mesh.cells[idx].x = i * mesh.dx;
            mesh.cells[idx].y = j * mesh.dy;
 
            // Define crooked pipe: horizontal then vertical bend
            bool in_pipe = (j >= 8 && j <= 12 && i < 30) || (i >= 28 && i <= 32 && j >= 8 && j <= 18);
            mesh.cells[idx].in_pipe = in_pipe;
            mesh.cells[idx].material_id = in_pipe ? 1 : 0;
        }

    }

 

    // Define boundaries

    for (int i = 0; i < mesh.num_cells; ++i) {
        if (mesh.cells[i].in_pipe) {
            if (is_wall_cell(mesh, i)) {
                BoundaryCondition bc;
                bc.cell = i;
                bc.type = WALL;
                mesh.boundaries.push_back(bc);
            }
        }
    }

    return mesh;

}

 

bool is_wall_cell(const Mesh& mesh, int idx) {
    int i = idx % mesh.nx;
    int j = idx / mesh.nx;
 
    for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
            if (di == 0 && dj == 0) continue;
            int ni = i + di;
            int nj = j + dj;

            if (ni < 0 || ni >= mesh.nx || nj < 0 || nj >= mesh.ny) continue;
            int nidx = nj * mesh.nx + ni;
            if (!mesh.cells[nidx].in_pipe) return true;
        }
    }

    return false;

}