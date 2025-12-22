#include <gtest/gtest.h>

#include "../../../include/solver.h"
#include "../../../include/physics.h"
#include "../../../include/geometry.h"
#include "../../../include/setup.h"



TEST(RadiationIntegration, PipeTwoMaterials) {


    //MPI_Init(NULL, NULL);
    
    Pars pars= Pars();
    pars.tini=2000.0;
    pars.nx=200;
    pars.ny=20; 
    pars.dx=0.1;
    pars.dy=0.1;
    pars.nz=1;

    // Setup: mesh and state with top-hat initial condition
    Mesh mesh(200, 20,0.1,0.1); // example dimensions
    

    
     // Define two materials: optically thin vs thick
    //Materials thin(0.01, 1.0);   // low absorption
    //Materials thick(10.0, 1.0);  // high absorption
    Materials materials=initialize_materials(mesh);
    MaterialProperties copper{0.15, 385.0};
    int status = materials.add_material(2, copper);  //add absorbant copper to right hand side



    // Assign materials along pipe
    /*for (int i = 0; i < 100; ++i)
        mesh.setMaterial(i, thin);
    for (int i = 100; i < 200; ++i)
        mesh.setMaterial(i, thick);   */

    for(int i=0;i<mesh.num_cells;i++) {
        int ix=i%mesh.nx;
        int iy=i/mesh.nx;
        if((ix>=mesh.nx/2) && mesh.cells[i].in_pipe) 
            mesh.cells[i].material_id=2; //copper  
    }



    State state = initialize_physics(mesh, materials,pars);  //stores the initial step
    // Initialize radiation source at left boundary
    state.getTemperatureField()[10] = 500.0;
    double tini=pars.tini;

    for(int icell=0; icell<mesh.num_cells; icell++) {
        int i= icell % pars.nx;
        int j= (icell / pars.nx) % pars.ny; 
        int k= icell /(pars.nx*pars.ny);

        if(i<pars.nx/2)
            //tini=pars.tmin+pars.tini*std::exp(-((i-2)*(i-2)+(j-pars.ny/2)*(j-pars.ny/2))/(4.0)); 
            tini=pars.tini;    
        else
           tini=300; //very cold?
        
        state.temperature[icell]=tini;   
    }



    //intialize solver for each state
    pars.dt=1e-6; //small time step for stability
    pars.emisscale=10.00;
    RadSolve solver(mesh.nx, mesh.ny,pars.nz, pars);
    //initialise the next radiation flux to the current flux 
    int i;    
    for(i=0; i<1000; i++)
    {
        
        if( (i%100)==0 )
            state.write_vtk_file(state.temperature, i, pars);
            //state.write_vtk_file(state.radiation_fluxn[0], i, pars);

        solver.solveRadiationTransport(mesh, state, pars, pars.time);   //FIXME
        solver.solveRadiationTransport(mesh, state, pars, pars.time);   //FIXME
        std::cout<<"solved radtrans"<<std::endl; 
        ;//solver.apply_milne_boundary_conditions(mesh, state, pars);   //CHECKME
        std::cout<<"applied milne bc"<<std::endl;
        solve_material_heating(mesh, state,pars);   //FIXME  // generates nans! in the copper region
        solver.UpdateBEmission(mesh, state, pars);  
        solver.UpdateRadFlux(mesh, state, pars);
        //solver.updatestate(pars,mesh,state);
        pars.time+=pars.dt;
    }
    

    // Check: radiation penetrates further in thin region than thick
    double thin_region_val  = state.getTemperatureField()[50+10*200];
    double thick_region_val = state.getTemperatureField()[150+10*200];

    std::cout<<"Thin region temp: "<<thin_region_val<<std::endl;
    std::cout<<"Thin region temp: "<<thick_region_val<<std::endl;
    MPI_Finalize();
    //EXPECT_GT(thin_region_val, thick_region_val);
    EXPECT_EQ(2, i);
}