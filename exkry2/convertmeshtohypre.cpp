#include <iostream>
#include <fstream>
#include <vector>
#include "mpi.h"
#include "HYPRE.h"
#include "HYPRE_IJ_mv.h"

// Define mesh storage
std::vector<std::vector<int>> elements;
std::vector<std::vector<double>> nodes;

// Read Gmsh .msh file (triangular elements)
void readMesh(const std::string& filename) {
    std::ifstream meshFile(filename);
    std::string line;
    bool readingNodes = false, readingElements = false;

    while (std::getline(meshFile, line)) {
        if (line.find("$Nodes") != std::string::npos) {
            readingNodes = true;
            continue;
        }
        if (line.find("$EndNodes") != std::string::npos) {
            readingNodes = false;
            continue;
        }
        if (line.find("$Elements") != std::string::npos) {
            readingElements = true;
            continue;
        }
        if (line.find("$EndElements") != std::string::npos) {
            readingElements = false;
            continue;
        }

        if (readingNodes) {
            int id;
            double x, y, z;
            std::stringstream ss(line);
            ss >> id >> x >> y >> z;
            nodes.push_back({x, y});
        }

        if (readingElements) {
            int id, type, tag, node1, node2, node3;
            std::stringstream ss(line);
            ss >> id >> type >> tag >> node1 >> node2 >> node3;
            if (type == 2) { // Type 2 is a triangle element in Gmsh
                elements.push_back({node1-1, node2-1, node3-1});
            }
        }
    }
    meshFile.close();
}

// Build HYPRE IJ Matrix from triangular mesh
void buildHypreMatrix(HYPRE_IJMatrix &A, int numNodes) {
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, numNodes-1, 0, numNodes-1, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    for (const auto& elem : elements) {
        int row = elem[0];
        int cols[3] = {elem[0], elem[1], elem[2]};
        double values[3] = {1.0, -0.5, -0.5};
        HYPRE_IJMatrixSetValues(A, 1, 3, &row, cols, values);
    }

    HYPRE_IJMatrixAssemble(A);
}

// Solve using Hypre PCG solver
void solveHypreSystem(HYPRE_IJMatrix &A, int numNodes) {
    HYPRE_IJVector b, x;
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, numNodes-1, &b);
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, numNodes-1, &x);

    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);
    HYPRE_IJVectorInitialize(x);

    for (int i = 0; i < numNodes; i++) {
        double rhs = (i == 0) ? 1.0 : 0.0; // Source at one node
        HYPRE_IJVectorSetValues(b, 1, &i, &rhs);
    }

    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorAssemble(x);

    // Setup and solve PCG
    HYPRE_Solver solver;
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b, par_x;
    
    HYPRE_IJMatrixGetObject(A, (void**) &par_A);
    HYPRE_IJVectorGetObject(b, (void**) &par_b);
    HYPRE_IJVectorGetObject(x, (void**) &par_x);

    HYPRE_PCGCreate(&solver);
    HYPRE_PCGSetup(solver, par_A, par_b, par_x);
    HYPRE_PCGSolve(solver, par_A, par_b, par_x);

    HYPRE_PCGDestroy(solver);
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " mesh.msh\n";
        MPI_Finalize();
        return 1;
    }

    readMesh(argv[1]);
    int numNodes = nodes.size();
    HYPRE_IJMatrix A;
    buildHypreMatrix(A, numNodes);
    solveHypreSystem(A, numNodes);

    MPI_Finalize();
    return 0;
}