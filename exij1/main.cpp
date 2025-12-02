#include <cstdio>
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_IJ_mv.h"

int main(int argc, char *argv[])
{
    // Initialize MPI
    int num_procs, myid;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // Problem size
    int n = 2; // 2x2 system
    int ilower = myid * n;
    int iupper = (myid + 1) * n - 1;

    // Create the matrix A
    HYPRE_IJMatrix A;
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);

    // Fill matrix entries (simple 2x2 system: [[4,1],[1,3]])
    int rows[2] = {ilower, ilower+1};
    int ncols[2] = {2,2};
    int cols[2][2] = {{ilower, ilower+1}, {ilower, ilower+1}};
    double values[2][2] = {{4.0, 1.0}, {1.0, 3.0}};

    for (int i = 0; i < 2; i++) {
        HYPRE_IJMatrixSetValues(A, 1, &ncols[i], &rows[i], cols[i], values[i]);
    }

    HYPRE_IJMatrixAssemble(A);

    // Create vectors b and x
    HYPRE_IJVector b;
    HYPRE_IJVector x;
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);
    HYPRE_IJVectorInitialize(x);

    // RHS b = [1,2]
    double rhs[2] = {1.0, 2.0};
    HYPRE_IJVectorSetValues(b, 2, rows, rhs);

    // Initial guess x = [0,0]
    double x0[2] = {0.0, 0.0};
    HYPRE_IJVectorSetValues(x, 2, rows, x0);

    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorAssemble(x);

    // Extract ParCSR objects
    HYPRE_ParCSRMatrix parA;
    HYPRE_ParVector parb, parx;
    HYPRE_IJMatrixGetObject(A, (void**) &parA);
    HYPRE_IJVectorGetObject(b, (void**) &parb);
    HYPRE_IJVectorGetObject(x, (void**) &parx);

    // Create solver (BoomerAMG)
    HYPRE_Solver solver;
    HYPRE_BoomerAMGCreate(&solver);
    HYPRE_BoomerAMGSetPrintLevel(solver, 2); // print residuals
    HYPRE_BoomerAMGSetMaxIter(solver, 20);

    // IMPORTANT: Setup phase
    HYPRE_BoomerAMGSetup(solver, parA, parb, parx);



    // Solve Ax = b
    HYPRE_BoomerAMGSolve(solver, parA, parb, parx);

    // Get solution
    double sol[2];
    HYPRE_IJVectorGetValues(x, 2, rows, sol);

    if (myid == 0) {
        printf("Solution: x = [%f, %f]\n", sol[0], sol[1]);
    }

    // Cleanup
    HYPRE_BoomerAMGDestroy(solver);
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);

    MPI_Finalize();
    return 0;
}