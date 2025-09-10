
#include <mpi.h>
#include "../include/geometry.h"
#include "../include/material.h"
#include "../include/physics.h"
#include "../include/solver.h"
#include "../include/boundary.h"

 

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

 

    // Setup

    Mesh mesh = setup_crooked_pipe_geometry();  //defined in geometry.h
    Materials materials = initialize_materials(mesh);

    State state = initialize_physics(mesh, materials);  //stores the initial step
    State state1(state);
    State state2(state);
    State statef(state); //1st half step update
 

    for (int timestep = 0; timestep < NSTEP; ++timestep) {

        apply_milne_boundary_conditions(mesh, state);

        linearize_emissive_source(mesh,state);  //see physics.h

        solve_radiation_groups(mesh, state);

        solve_material_heating(mesh, state);

    }

 

    MPI_Finalize();

    return 0;

}