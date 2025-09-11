#include "../include/solver.h"
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"

 


int setupsolver() {
    // This function sets up the solver parameters and initializes necessary variables
    // It can be expanded to include more complex setup logic if needed
    // For now, it simply returns true to indicate successful setup
    return 0;
}

void solve_radiation_groups(const Mesh& mesh, State& state) {

    for (int g = 0; g < NUM_GROUPS; ++g) {

        // Setup HYPRE solver for group g

        HYPRE_StructGrid grid;
        HYPRE_StructMatrix A;
        HYPRE_StructVector b, x;
        HYPRE_StructSolver solver;

 

        // Initialize grid, matrix, vectors (omitted for brevity)

 

        // Fill matrix A and vector b using diffusion equation

        // A[i][j] = -D_g * Laplacian + sigma_a

        // b[i] = source_term[g][i]

 

        // Solve

        HYPRE_StructGMRESCreate(MPI_COMM_WORLD, &solver);
        HYPRE_StructGMRESSetTol(solver, 1e-6);
        HYPRE_StructGMRESSetup(solver, A, b, x);
        HYPRE_StructGMRESSolve(solver, A, b, x);

 

        // Extract solution into state.radiation_flux[g]

        // (omitted for brevity)

 

        HYPRE_StructGMRESDestroy(solver);

    }

}


    //RadSolve::RadSolve(int mnx, int mny) : nx(mnx), ny(mny), T(mnx, std::vector<double>(mny, 300.0)) {
    RadSolve::RadSolve(int mnx, int mny) : nx(mnx), ny(mny) {
        MPI_Init(NULL, NULL);
        HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
        HYPRE_StructStencilCreate(2, 5, &stencil);
        HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
        HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
    }

    RadSolve::~RadSolve() {
        HYPRE_StructMatrixDestroy(A);
        HYPRE_StructGridDestroy(grid);
        HYPRE_StructStencilDestroy(stencil);
        HYPRE_StructVectorDestroy(b);
        HYPRE_StructVectorDestroy(x);
        MPI_Finalize();
    }

    // Method to return the temperature field
    /*const std::vector<std::vector<double>>& RadSolve::getTemperatureField() const {
        return T;
    }*/

    // Example function to modify the temperature field
    /*void RadSolve::setTemperature(int i, int j, double value) {
        if (i >= 0 && i < nx && j >= 0 && j < ny) {
            T[i][j] = value;
        } else {
            std::cerr << "Error: Index out of bounds!" << std::endl;
        }
    }*/



    /*void RadSolve::readMesh(const std::string& filename) {
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
    }*/

    /*void RadSolve::setupGrid() {
        int ilower[2] = {0, 0};
        int iupper[2] = {nx-1, ny-1};
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
    }*/

    void RadSolve::solveRadiationTransport(const Mesh& mesh, State& state, double t) {
        //double T[NX][NY];
        /*for (int i = 0; i < nx; i++)
            for (int j = 0; j < ny; j++)
                T[i][j] = 300.0;  // Initial temperature field*/

        /*for (int t = 0; t < TIME_STEPS; t++) {*/
            std::cout << "Solving time step " << t << std::endl;

            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
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

            for (int i = 0; i < nx; i++) {
                for (int j = 0; j < ny; j++) {
                    int index[2] = {i, j};
                    double sol_value;
                    HYPRE_StructVectorGetValues(x, index, &sol_value);
                    T[i][j] = sol_value;
                }
            }

            HYPRE_StructPCGDestroy(solver);
        /*}*/
    }
