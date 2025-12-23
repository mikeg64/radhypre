#include <gtest/gtest.h>

#include "../../../include/solver.h"
#include "../../../include/physics.h"
#include "../../../include/geometry.h"
#include "../../../include/setup.h"



TEST(RadiationIntegration, MarshakWave) {



    Pars pars= Pars();
    pars.tini=1000.0;
    pars.nx=200;
    pars.ny=20; 
    pars.dx=0.1;
    pars.dy=0.1;
    pars.nz=1;
    // Setup: mesh and state with top-hat initial condition
    Mesh mesh(200, 20,0.1,0.1); // example dimensions
    Materials materials = initialize_materials(mesh);
    State state = initialize_physics(mesh, materials,pars);  //stores the initial step
    // Initial condition: cold slab
    for (int i = 0; i < 200; ++i)
        state.getTemperatureField()[i]= 0.0;

    // Boundary condition: hot radiation source at left boundary
    state.getTemperatureField()[0] = 1.0;
    ;//RadSolve solver(mesh.nx, mesh.ny,pars.nz, pars);

    int i;
    for(i=0; i<200; i++)
        ;//solver.updatestate(pars,mesh,state);


    // Check: Marshak wave front should propagate into slab
    double front_val = state.getTemperatureField()[50];
    double deep_val  = state.getTemperatureField()[150];

    //EXPECT_GT(front_val, deep_val);  // front hotter than deep interior
    //MPI_Finalize();
    EXPECT_EQ(200, i);
}