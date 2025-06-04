#include <Kokkos_Core.hpp>

#include "HYPRE.h"

#include "HYPRE_parcsr_ls.h"

 

int main(int argc, char *argv[]) {

    Kokkos::initialize(argc, argv);

    {

        const int nx = 128, ny = 128;

 

        Kokkos::View<double**> u("solution", nx, ny);

        Kokkos::View<double**> rhs("rhs", nx, ny);

 

        Kokkos::parallel_for("init_rhs", nx * ny, KOKKOS_LAMBDA(int idx) {

            int i = idx % nx;

            int j = idx / nx;

            rhs(i, j) = 1.0; // Example source term

            u(i, j) = 0.0;

        });

 

        HYPRE_IJMatrix A;

        HYPRE_IJVector b, x;

        // Initialize and fill A, b, x...

 

        HYPRE_ParCSRMatrix par_A;

        HYPRE_ParVector par_b, par_x;

        HYPRE_IJMatrixGetObject(A, (void**) &par_A);

        HYPRE_IJVectorGetObject(b, (void**) &par_b);

        HYPRE_IJVectorGetObject(x, (void**) &par_x);

 

        HYPRE_Solver solver;

        HYPRE_BoomerAMGCreate(&solver);

        HYPRE_BoomerAMGSetPrintLevel(solver, 1);

        HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x);

        HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x);

 

        HYPRE_BoomerAMGDestroy(solver);

        HYPRE_IJMatrixDestroy(A);

        HYPRE_IJVectorDestroy(b);

        HYPRE_IJVectorDestroy(x);

    }

    Kokkos::finalize();

    return 0;

}