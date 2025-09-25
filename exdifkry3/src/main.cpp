
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
    Mesh mesh = setup_crooked_pipe_geometry(pars);  //defined in geometry.h
    Materials materials = initialize_materials(mesh);

    State state = initialize_physics(mesh, materials);  //stores the initial step
    State state1(state);
    State statef(state);
    State state2(state);
    

    //intialize solver for each state
    RadSolve solver(mesh.nx, mesh.ny,pars);
 
    double temptol=1e-3; //temperature tolerance for convergence
    double maxtempdif=0.0;
    double time=0.0;
    double dt=1e-10; //initial time step
    // Main time-stepping loop
    for (int timestep = 0; timestep < pars.nstep; ++timestep) {
        state1.copy(state);
        updatestate(pars, mesh, state, state1, state2, statef, solver);
        

        pars.time+=dt;
        std::cout<<"Completed timestep "<<timestep<<" time "<<time<<" dt "<<dt<< std::endl;
        //maxtempdif<0.5*temptol then new dt=2*dt update state to state2

        //if maxtempdif>temptol then new dt=0.5*dt update state to state1
        if(timestep%pars.nsaveinterval==0)
            state.write_vtk_file(state.temperature, (1+(timestep/pars.nsaveinterval)), pars);

    }

 

    MPI_Finalize();

    return 0;

}

int updatestate(Pars &pars, Mesh &mesh, State &state, State &state1, State &state2, State &statef, RadSolve &solver)
{
    int status=0;    
    
    //take a full step dt
        
        
        //solve_radiation_groups(mesh, state);
        solver.solveRadiationTransport(mesh, state, pars, pars.time);
        //linearize_emissive_source(mesh,state,pars);  //see physics.h   
        solver.apply_milne_boundary_conditions(mesh, state, pars);
        solve_material_heating(mesh, state,pars);

        //Take 2 half steps dt=dt/2

        pars.dt=0.5*pars.dt;
        solver.solveRadiationTransport(mesh, state1, pars, pars.time);
        //solve_radiation_groups(mesh, state1);   
        solver.apply_milne_boundary_conditions(mesh, state1, pars);
        solve_material_heating(mesh, state1,pars);
        //linearize_emissive_source(mesh,state1,pars);  //see physics.h

        state2.copy(state1);
        solver.solveRadiationTransport(mesh, state2, pars, pars.time);
        //solve_radiation_groups(mesh, state2);
        solver.apply_milne_boundary_conditions(mesh, state2, pars);
        solve_material_heating(mesh, state2,pars);
        //linearize_emissive_source(mesh,state2,pars);  //see physics.h
        pars.dt=2.0*pars.dt;

        //put this in a routine called convergence check and update dt
        //compare state and state2
        //eg something like this
        double maxdifrat=0.0;
        double dif,difrat;
        for(int i=0;i<mesh.num_cells;i++) {
            dif=fabs(state.temperature[i]-state2.temperature[i]);
            difrat=dif/(0.5*(state.temperature[i]+state2.temperature[i])+1e-10);
            maxdifrat=(difrat>maxdifrat?difrat:maxdifrat);
        }
        //do some trickery to reduce the timestep
        if(maxdifrat<pars.temptol && pars.dt<pars.dtmax) {
            pars.dt=2.0*pars.dt;
            state.copy(state2);
        }
        else if(maxdifrat>2.0*pars.temptol && pars.dt>pars.dtmin)
         {
            pars.dt=0.5*pars.dt;
            state.copy(state1);
        }
        //else keep time step the same and just use state

        return status;
}
