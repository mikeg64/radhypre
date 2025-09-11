
#include <mpi.h>
#include "../include/geometry.h"
#include "../include/material.h"
#include "../include/physics.h"
#include "../include/solver.h"
#include "../include/boundary.h"
#include "../include/setup.h"

int updatestate(Pars &pars, Mesh &mesh, State &state, State &state1, State &state2, State &statef, RadSolve &solver); 

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    HYPRE_Init();

 

    // Setup
    Pars pars= Pars();
    Mesh mesh = setup_crooked_pipe_geometry();  //defined in geometry.h
    Materials materials = initialize_materials(mesh);

    State state = initialize_physics(mesh, materials);  //stores the initial step
    State state1(state);
    State statef(state);
    State state2(state);
    

    //intialize solver for each state
    RadSolve solver(mesh.nx, mesh.ny);
 
    double temptol=1e-3; //temperature tolerance for convergence
    double maxtempdif=0.0;
    double time=0.0;
    double dt=1e-10; //initial time step
    // Main time-stepping loop
    for (int timestep = 0; timestep < NSTEP; ++timestep) {
        state1.copy(state);
        updatestate(pars, mesh, state, state1, state2, statef, solver);
        

        time+=dt;
        std::cout<<"Completed timestep "<<timestep<<" time "<<time<<" dt "<<dt<< std::endl;
        //maxtempdif<0.5*temptol then new dt=2*dt update state to state2

        //if maxtempdif>temptol then new dt=0.5*dt update state to state1


    }

 

    MPI_Finalize();

    return 0;

}

int updatestate(Pars &pars, Mesh &mesh, State &state, State &state1, State &state2, State &statef, RadSolve &solver)
{
        //take a full step dt
        apply_milne_boundary_conditions(mesh, state);
        linearize_emissive_source(mesh,state);  //see physics.h
        solve_radiation_groups(mesh, state);
        solver.solveRadiationTransport(mesh, state, pars.time);
        solve_material_heating(mesh, state);

        //Take 2 half steps dt=dt/2
        apply_milne_boundary_conditions(mesh, state1);
        linearize_emissive_source(mesh,state1);  //see physics.h
        solver.solveRadiationTransport(mesh, state, pars.time);
        solve_radiation_groups(mesh, state1);
        solve_material_heating(mesh, state1);

        state2.copy(state1);
        apply_milne_boundary_conditions(mesh, state2);
        linearize_emissive_source(mesh,state2);  //see physics.h
        solver.solveRadiationTransport(mesh, state, pars.time);
        solve_radiation_groups(mesh, state2);
        solve_material_heating(mesh, state2);

        //put this in a routine called convergence check and update dt
        //compare state and state2
        //eg something like this
        for(int i=0;i<mesh.num_cells;i++) {
            if(fabs(state.temperature[i]-state2.temperature[i])>pars.temptol) {
                std::cout<<"timestep "<<pars.dt<<" cell "<<i<<" temp1 "<<state.temperature[i]<<" temp2 "<<state2.temperature[i]<<std::endl;
            }
            statef.temperature[i]=0.5*(state.temperature[i]+state2.temperature[i]);
        }
        //do some trickery to reduce the timestep
}
