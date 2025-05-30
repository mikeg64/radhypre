
#include <iostream>
#include <cmath>
#include <fstream>
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_krylov.h"
#include "mpi.h"

#define NX 100
#define NY 100
#define NZ 100
#define L 1.0

// Function to calculate diffusion coefficient D(x, y, z)
double D(double x, double y, double z) {
    return 0.1 + 0.05 * sin(2 * M_PI * x) * sin(2 * M_PI * y) * sin(2 * M_PI * z);
}

// Function to calculate absorption coefficient kappa(x, y, z)
double kappa(double x, double y, double z) {
    return 1.0 + 0.5 * cos(2 * M_PI * x) * cos(2 * M_PI * y) * cos(2 * M_PI * z);
}

// Function to calculate source term S(x, y, z)
double S(double x, double y, double z) {
    double sigma = 0.01;
    double x0 = 0.5, y0 = 0.5, z0 = 0.1;
    return exp(-((x - x0) * (x - x0) + (y - y0) * (y - y0) + (z - z0) * (z - z0)) / (2 * sigma * sigma));
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int num_procs, my_id;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    int n = NX * NY * NZ;
    int ilower = my_id * n / num_procs;
    int iupper = (my_id + 1) * n / num_procs - 1;

    HYPRE_IJMatrix A;
    HYPRE_IJVector b;
    HYPRE_IJVector x;

    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    int *rows = new int[n];
    double *values = new double[n];
    double *rhs_values = new double[n];

    for (int i = 0; i < n; i++) {
        rows[i] = ilower + i;
        values[i] = 0.0;
        rhs_values[i] = 0.0;
    }

    double hx = L / (NX - 1);
    double hy = L / (NY - 1);
    double hz = L / (NZ - 1);

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int idx = i * NY * NZ + j * NZ + k;
                double x = i * hx;
                double y = j * hy;
                double z = k * hz;

                double D_val = D(x, y, z);
                double kappa_val = kappa(x, y, z);
                double S_val = S(x, y, z);

                values[idx] = -2.0 * D_val / (hx * hx) - 2.0 * D_val / (hy * hy) - 2.0 * D_val / (hz * hz) + kappa_val;
                rhs_values[idx] = S_val;

                if (i > 0) {
                    int idx_left = (i - 1) * NY * NZ + j * NZ + k;
                    values[idx_left] = D_val / (hx * hx);
                }
                if (i < NX - 1) {
                    int idx_right = (i + 1) * NY * NZ + j * NZ + k;
                    values[idx_right] = D_val / (hx * hx);
                }
                if (j > 0) {
                    int idx_down = i * NY * NZ + (j - 1) * NZ + k;
                    values[idx_down] = D_val / (hy * hy);
                }
                if (j < NY - 1) {
                    int idx_up = i * NY * NZ + (j + 1) * NZ + k;
                    values[idx_up] = D_val / (hy * hy);
                }
                if (k > 0) {
                    int idx_back = i * NY * NZ + j * NZ + (k - 1);
                    values[idx_back] = D_val / (hz * hz);
                }
                if (k < NZ - 1) {
                    int idx_front = i * NY * NZ + j * NZ + (k + 1);
                    values[idx_front] = D_val / (hz * hz);
                }
            }
        }
    }

    HYPRE_IJMatrixSetValues(A, n, NULL, rows, values);
    HYPRE_IJVectorSetValues(b, n, rows, rhs_values);

    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorAssemble(x);

    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;

    HYPRE_IJMatrixGetObject(A, (void **) &par_A);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);

    HYPRE_Solver solver;
    HYPRE_BoomerAMGCreate(&solver);
    HYPRE_BoomerAMGSetPrintLevel(solver, 2);
    HYPRE_BoomerAMGSetMaxIter(solver, 100);
    HYPRE_BoomerAMGSetTol(solver, 1e-6);

    HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x);
    HYPRE_BoomerAMGIterate(solver, par_A, par_b, par_x);

    double *solution = new double[n];
    HYPRE_IJVectorGetValues(x, n, rows, solution);

    std::ofstream outfile("solution.txt");
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                int idx = i * NY * NZ + j * NZ + k;
                outfile << i << " " << j << " " << k << " " << solution[idx] << std::endl;
            }
        }
    }
    outfile.close();

    delete[] rows;
    delete[] values;
    delete[] rhs_values;
    delete[] solution;

    HYPRE_BoomerAMGDestroy(solver);
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);

    MPI_Finalize();
    return 0;
}
