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
                bool in_pipe =(j>=1 && j<=(ny-2) && i>=1 && i<=(nx-2)); //temp for box test           
                cells[idx].in_pipe = in_pipe;
                cells[idx].material_id = (in_pipe ? 1 : 0);   // the first material 0  is the wall material absorber/reflector high heat capacity
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




Mesh setup_crooked_pipe_geometry(Pars &pars) {
    Mesh mesh(pars.nx, pars.ny, pars.dx, pars.dy);
    
    //mesh.num_cells = mesh.nx * mesh.ny;
    //mesh.cells.resize(mesh.num_cells);

 
    int ibound=0;
    for (int j = 0; j < mesh.ny; ++j) {
        for (int i = 0; i < mesh.nx; ++i) {
            int idx = j * mesh.nx + i;
            mesh.cells[idx].x = i * mesh.dx;
            mesh.cells[idx].y = j * mesh.dy;
 
            // Define crooked pipe: horizontal then vertical bend
           // bool in_pipe = (j >= 8 && j <= 12 && i < 30) || (i >= 28 && i <= 32 && j >= 8 && j <= 18);
            bool in_pipe = (j >= (NY/5) && j <= (2*NY/5) && i < (2*NX/7)) 
                               || (i >= (2*NX/7) && i <= (3*NX/7) && j >= (NY/5) && j <= (4*NY/5))
                               || (j >= (3*NY/5) && j <= (4*NY/5) && i >= (3*NX/7) && i <= (4*NX/7)) 
                               || (i >= (4*NX/7) && i <= (5*NX/7) && j >= (NY/5) && j <= (4*NY/5))
                               || (j >= (NY/5) && j <= (2*NY/5) && i >= (5*NX/7) && i <= (NX));




                               
            mesh.cells[idx].in_pipe = in_pipe;
            mesh.cells[idx].material_id = in_pipe ? 1 : 0;

            if(in_pipe)
            {
                if (is_wall_cell(mesh, idx)) {
                    BoundaryCondition bc;
                    bc.cell = idx;
                    bc.type = WALL;
                    // Determine wall type based on position
                    // Left wall, right wall, upper wall, lower wall

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

                    bc.wall_type = (WallType)ibound;
                    mesh.boundaries.push_back(bc);
                }

            }
            else {
                if(mesh.cells[idx].x<(2*mesh.dx) && mesh.cells[idx].y<(2*mesh.dy*mesh.ny/5)   && mesh.cells[idx].y> (mesh.dy*mesh.ny/5)  ) 
                { // Top row
                    BoundaryCondition bc;
                    bc.cell = idx;
                    bc.type = OUTLET;
                    mesh.boundaries.push_back(bc);
                }
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