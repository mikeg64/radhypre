#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "mpi.h"
#include "HYPRE.h"
#include "HYPRE_struct_ls.h"
#include "ex.h"

#define NX 50  // Grid points in X
#define NY 50  // Grid points in Y
#define TIME_STEPS 10 // Number of time steps
#define DT 0.01 // Time step size
#define ALPHA 0.1 // Absorption coefficient

class RadSolve {
public:
    RadSolve() {
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
        double T[NX][NY];
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

private:
    HYPRE_StructGrid grid;
    HYPRE_StructStencil stencil;
    HYPRE_StructMatrix A;
    HYPRE_StructVector b, x;
    HYPRE_StructSolver solver;
    std::vector<std::pair<double, double>> nodes;
};

int main(int argc, char *argv[]) {
    RadSolve solver;
    solver.readMesh("mesh.msh");  // Read mesh from Gmsh file
    solver.setupGrid();
    solver.solveRadiationTransport();
    return 0;
}