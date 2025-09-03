
#include <mpi.h>
#include "../include/geometry.h"
#include "materials.h"
#include "physics.h"
#include "solver.h"
#include "boundary.h"

 

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

 

    // Setup

    Mesh mesh = setup_crooked_pipe_geometry();
    Materials materials = initialize_materials(mesh);
    PhysicsState state = initialize_physics(mesh, materials);

 

    for (int timestep = 0; timestep < MAX_TIMESTEPS; ++timestep) {

        apply_milne_boundary_conditions(mesh, state);

        linearize_emissive_source(state);

        solve_radiation_groups(mesh, state);

        solve_material_heating(mesh, state);

    }

 

    MPI_Finalize();

    return 0;

}