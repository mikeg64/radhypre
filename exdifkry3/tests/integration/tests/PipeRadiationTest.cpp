#include <gtest/gtest.h>

#include "../../../include/solver.h"
#include "../../../include/physics.h"
#include "../../../include/geometry.h"
#include "../../../include/setup.h"



TEST(RadiationIntegration, PipeTwoMaterials) {


    Pars pars= Pars();
    // Setup: mesh and state with top-hat initial condition
    Mesh mesh(200, 20,0.1,0.1); // example dimensions
    //Materials materials = initialize_materials(mesh);
   
    Materials materials=initialize_materials(mesh);
    // Define two materials: optically thin vs thick
    //Materials thin(0.01, 1.0);   // low absorption
    //Materials thick(10.0, 1.0);  // high absorption

    // Assign materials along pipe
    /*for (int i = 0; i < 100; ++i)
        mesh.setMaterial(i, thin);
    for (int i = 100; i < 200; ++i)
        mesh.setMaterial(i, thick);   */

     State state = initialize_physics(mesh, materials,pars);  //stores the initial step
    // Initialize radiation source at left boundary
    state.getTemperatureField()[10] = 1.0;
    //intialize solver for each state
    RadSolve solver(mesh.nx, mesh.ny,pars.nz, pars);

    for(int i=0; i<100; i++)
        solver.updatestate(pars,mesh,state);

    

    // Check: radiation penetrates further in thin region than thick
    double thin_region_val  = state.getTemperatureField()[50+10*200];
    double thick_region_val = state.getTemperatureField()[150+10*200];

    EXPECT_GT(thin_region_val, thick_region_val);
}