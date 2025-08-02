#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_IJ_mv.h>
#include <cmath>
#include <vector>
#include <iostream>

const int N = 100; // Number of grid points
const double L = 1.0; // Domain length
const double dx = L / (N - 1);
const double kappa0 = 10.0; // Opacity inside top-hat
const double source = 1.0; // Source function inside top-hat

// Top-hat region
const double x_min = 0.4;
const double x_max = 0.6;

double kappa(double x) {
    return (x >= x_min && x <= x_max) ? kappa0 : 0.0;
}

double S(double x) {
    return (x >= x_min && x <= x_max) ? source : 0.0;
}

int main() {
    HYPRE_IJMatrix A;
    HYPRE_IJVector b, x;
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b, par_x;
    HYPRE_Solver solver;

    // Initialize HYPRE
    HYPRE_Init();

    // Create matrix and vectors
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N-1, 0, N-1, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    // Fill matrix and RHS
    for (int i = 0; i < N; ++i) {
        double xi = i * dx;
        double k = kappa(xi);
        double s = S(xi);

        int row = i;
        int ncols = 1;
        double val = 1.0 + k * dx;
        int col = i;
        HYPRE_IJMatrixSetValues(A, 1, &ncols, &row, &col, &val);

        double rhs = k * s * dx;
        HYPRE_IJVectorSetValues(b, 1, &row, &rhs);
        double x0 = 0.0;
        HYPRE_IJVectorSetValues(x, 1, &row, &x0);
    }

    // Finalize
    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJMatrixGetObject(A, (void**) &par_A);

    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void**) &par_b);

    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void**) &par_x);

    // Create and setup PCG solver
    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_PCGSetMaxIter(solver, 100);
    HYPRE_PCGSetTol(solver, 1e-8);
    HYPRE_ParCSRPCGSetup(solver, par_A, par_b, par_x);
    HYPRE_ParCSRPCGSolve(solver, par_A, par_b, par_x);

    // Cleanup
    HYPRE_ParCSRPCGDestroy(solver);
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);
    HYPRE_Finalize();

    return 0;
}