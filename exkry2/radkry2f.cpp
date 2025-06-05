#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "mpi.h"
#include "HYPRE.h"
#include "HYPRE_struct_ls.h"
#include "ex.h"

#define GNX 50  // Grid points in X
#define GNY 50  // Grid points in Y
#define TIME_STEPS 10 // Number of time steps
#define DT 0.01 // Time step size
#define ALPHA 0.1 // Absorption coefficient

class RadSolve {
public:
    RadSolve(int nx, int ny) : NX(nx), NY(ny), T(nx, std::vector<double>(ny, 300.0)) {
        MPI_Init(NULL, NULL);
        HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
        HYPRE_StructStencilCreate(2, 5, &stencil);
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    }

    ~RadSolve() {
        HYPRE_StructMatrixDestroy(A);
        HYPRE_StructGridDestroy(grid);
        HYPRE_StructStencilDestroy(stencil);
        HYPRE_StructVectorDestroy(b);
        HYPRE_StructVectorDestroy(x);
        MPI_Finalize();
    }

    // Method to return the temperature field
    const std::vector<std::vector<double>>& getTemperatureField() const {
        return T;
    }

    // Example function to modify the temperature field
    void setTemperature(int i, int j, double value) {
        if (i >= 0 && i < NX && j >= 0 && j < NY) {
            T[i][j] = value;
        } else {
            std::cerr << "Error: Index out of bounds!" << std::endl;
        }
    }



    void readMesh(const std::string& filename) {
        std::ifstream meshFile(filename);
        std::string line;
        while (std::getline(meshFile, line)) {
            if (line.find("Nodes") != std::string::npos) {
                while (std::getline(meshFile, line)) {
                    if (line.find("EndNodes") != std::string::npos) break;
                    double x, y;
                    int index;
                    sscanf(line.c_str(), "%d %lf %lf", &index, &x, &y);
                    nodes.push_back({x, y});
                }
            }
        }
        meshFile.close();
    }

    void setupGrid() {
        int ilower[2] = {0, 0};
        int iupper[2] = {NX-1, NY-1};
        HYPRE_StructGridSetExtents(grid, ilower, iupper);
        HYPRE_StructGridAssemble(grid);

        int offsets[5][2] = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}};
        for (int i = 0; i < 5; i++)
            HYPRE_StructStencilSetElement(stencil, i, offsets[i]);

           // Create vectors
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
        HYPRE_StructVectorInitialize(b);
        HYPRE_StructVectorInitialize(x);

        //HYPRE_StructMatrixSetStencil(A, stencil);
        HYPRE_StructMatrixInitialize(A);
    }

    void solveRadiationTransport() {
        //double T[NX][NY];
        for (int i = 0; i < NX; i++)
            for (int j = 0; j < NY; j++)
                T[i][j] = 300.0;  // Initial temperature field

        for (int t = 0; t < TIME_STEPS; t++) {
            std::cout << "Solving time step " << t << std::endl;

            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int index[2] = {i, j};
                    double source = (i == 0) ? 6000.0 : 0.0;  // Source at left boundary
                    double rhs_value = T[i][j] + DT * source;
                    HYPRE_StructVectorSetValues(b, index, rhs_value);
                }
            }

            HYPRE_StructVectorAssemble(b);
            HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
            HYPRE_StructPCGSetTol(solver, 1e-6);
            HYPRE_StructPCGSetMaxIter(solver, 100);
            HYPRE_StructPCGSetup(solver, A, b, x);
            HYPRE_StructPCGSolve(solver, A, b, x);

            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    int index[2] = {i, j};
                    double sol_value;
                    HYPRE_StructVectorGetValues(x, index, &sol_value);
                    T[i][j] = sol_value;
                }
            }

            HYPRE_StructPCGDestroy(solver);
        }
    }

   // Function to write temperature data to a VTK file
    void write_vtk_file(const std::vector<std::vector<double>> &T, int time_step, const std::string &filename = "temperature") {
        std::string vtk_filename = filename + "_" + std::to_string(time_step) + ".vtk";
        std::ofstream vtk_file(vtk_filename);

        if (!vtk_file) {
            std::cerr << "Error: Unable to open file for writing!" << std::endl;
            return;
        }

        vtk_file << "# vtk DataFile Version 3.0\n";
        vtk_file << "Temperature field output\n";
        vtk_file << "ASCII\n";
        vtk_file << "DATASET STRUCTURED_POINTS\n";
        vtk_file << "DIMENSIONS " << NX << " " << NY << " 1\n";
        vtk_file << "ORIGIN 0 0 0\n";
        vtk_file << "SPACING 1.0 1.0 1.0\n";
        vtk_file << "POINT_DATA " << NX * NY << "\n";
        vtk_file << "SCALARS Temperature float\n";
        vtk_file << "LOOKUP_TABLE default\n";

        // Write the temperature values in row-major order
        for (int j = 0; j < NY; j++) {
            for (int i = 0; i < NX; i++) {
                vtk_file << T[i][j] << "\n";
            }
        }

        vtk_file.close();
        std::cout << "VTK file '" << vtk_filename << "' written successfully!" << std::endl;
    }




private:
    int NX, NY; // Grid dimensions 
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    std::vector<std::vector<double>> T; // Temperature field
    HYPRE_StructSolver solver;
    std::vector<std::pair<double, double>> nodes;
};

int main(int argc, char *argv[]) {
    RadSolve solver(GNX, GNY);  // Create RadSolve instance with grid dimensions
    solver.readMesh("mesh.msh");  // Read mesh from Gmsh file
    solver.setupGrid();
    solver.solveRadiationTransport();
    solver.write_vtk_file(solver.getTemperatureField(), 0);  // Write initial temperature field to VTK file
    return 0;
}