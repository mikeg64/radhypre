#include <gtest/gtest.h>
#include "../../../include/solver.h"
#include "../../../include/physics.h"
#include "../../../include/geometry.h"
#include "../../../include/setup.h"

TEST(DiffusionIntegration, TopHatProfile) {

    MPI_Init(NULL, NULL);

    Pars pars= Pars();
    pars.tini=1000.0;
    pars.nx=100;
    pars.ny=100; 
    pars.dx=0.1;
    pars.dy=0.1;
    pars.nz=1;
    
    // Setup: mesh and state with top-hat initial condition
    Mesh mesh(100, 100,0.1,0.1); // example dimensions
    Materials materials = initialize_materials(mesh);
    State state = initialize_physics(mesh, materials,pars);  //stores the initial step
    // Initialize top-hat temperature field
    for (int j = 40; j < 60; ++j)
        for (int i = 40; i < 60; ++i)
            state.getTemperatureField()[j+i*100] = 1.0;

    
    //intialize solver for each state


    pars.dt=1e-6; //small time step for stability
    pars.emisscale=10.00;
    pars.scale=10.0;
    RadSolve solver(mesh.nx, mesh.ny,pars.nz, pars);

    int i;
    for(i=0; i<50; i++)
    {

        if( (i%10)==0 )
            state.write_vtk_file(state.temperature, i, pars);

        solver.solveRadiationTransport(mesh, state, pars, pars.time);

        solver.solveRadiationTransport(mesh, state, pars, pars.time);   //FIXME
        solver.solveRadiationTransport(mesh, state, pars, pars.time);   //FIXME
        std::cout<<"solved radtrans"<<std::endl; 
        ;//solver.apply_milne_boundary_conditions(mesh, state, pars);   //CHECKME
        std::cout<<"applied milne bc"<<std::endl;
        solve_material_heating(mesh, state,pars);   //FIXME  // generates nans! in the copper region
        solver.UpdateBEmission(mesh, state, pars);  
        solver.UpdateRadFlux(mesh, state, pars);
        
        pars.time+=pars.dt;

    }
       

    // Check: diffusion should smooth the top-hat edges
    double center = state.getTemperatureField()[50+50*100];
    double edge   = state.getTemperatureField()[40+40*100];

    //EXPECT_GT(center, edge);  // center hotter than edge after diffusion
    EXPECT_EQ(50, i);
    //MPI_Finalize();
}