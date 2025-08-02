const int Nx = 50, Ny = 50;



#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_IJ_mv.h>
#include <vector>
#include <cmath>
#include <iostream>

const int Nx = 100; // Number of grid points
const int Nx = 50; // Number of grid points


const double Lx = 1.0, Ly = 1.0;
const double dx = Lx / (Nx - 1), dy = Ly / (Ny - 1);
const double kappa0 = 10.0;
const double source = 1.0;
const double x_min = 0.4, x_max = 0.6;
const double y_min = 0.4, y_max = 0.6;

int idx(int i, int j, int Nx) {
    return i + j * Nx;
}

double kappa(double x, double y) {
    return (x >= x_min && x <= x_max && y >= y_min && y <= y_max) ? kappa0 : 0.0;
}

double S(double x, double y) {
    return (x >= x_min && x <= x_max && y >= y_min && y <= y_max) ? source : 0.0;
}

int main() {
    HYPRE_Init();

    int N = Nx * Ny;
    HYPRE_IJMatrix A;
    HYPRE_IJVector b, x;
    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b, par_x;
    HYPRE_Solver solver;

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, 0, N-1, 0, N-1, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, 0, N-1, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int row = idx(i, j, Nx);
            double xi = i * dx;
            double yj = j * dy;
            double k = kappa(xi, yj);
            double s = S(xi, yj);

            std::vector<int> cols;
            std::vector<double> vals;

            // Central point
            cols.push_back(row);
            vals.push_back(1.0 + 2.0 * k * dx);

            // Left neighbor
            if (i > 0) {
                cols.push_back(idx(i-1, j, Nx));
                vals.push_back(-0.5);
            }

            // Right neighbor
            if (i < Nx-1) {
                cols.push_back(idx(i+1, j, Nx));
                vals.push_back(-0.5);
            }

            // Bottom neighbor
            if (j > 0) {
                cols.push_back(idx(i, j-1, Nx));
                vals.push_back(-0.5);
            }

            // Top neighbor
            if (j < Ny-1) {
                cols.push_back(idx(i, j+1, Nx));
                vals.push_back(-0.5);
            }

            HYPRE_IJMatrixSetValues(A, 1, (int[]){(int)cols.size()}, &row, cols.data(), vals.data());

            double rhs = k * s * dx * dy;
            HYPRE_IJVectorSetValues(b, 1, &row, &rhs);
            double x0 = 0.0;
            HYPRE_IJVectorSetValues(x, 1, &row, &x0);
        }
    }

    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJMatrixGetObject(A, (void**) &par_A);

    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void**) &par_b);

    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void**) &par_x);

    HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
    HYPRE_PCGSetMaxIter(solver, 100);
    HYPRE_PCGSetTol(solver, 1e-8);
    HYPRE_ParCSRPCGSetup(solver, par_A, par_b, par_x);
    HYPRE_ParCSRPCGSolve(solver, par_A, par_b, par_x);

    HYPRE_ParCSRPCGDestroy(solver);
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);
    HYPRE_Finalize();

    return 0;
}