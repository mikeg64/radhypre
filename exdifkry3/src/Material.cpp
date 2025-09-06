#include    "../include/material.h"
#include    "../include/geometry.h"
#include    "../include/setup.h"



Materials initialize_materials(const Mesh& mesh)
{
    Materials mats;
    int num_cells = mesh.num_cells;
    //mats.sigma_a.resize(num_cells);
    //mats.heat_capacity.resize(num_cells);

    for (int i = 0; i < num_cells; ++i) {
        int material_id = mesh.cells[i].material_id;
        // Example: Assign properties based on material_id
        if (material_id == 1) { // Material 1
            mats.set_sigma_a(i, 100.0); // Example value
            mats.set_heat_capacity(i,  10.0); // Example value
        } else { // Default material
            mats.set_sigma_a(i, 10.0); // Example value
            mats.set_heat_capacity(i,  1.0); // Example value
        }
    }

    return mats;
}


/*Materials initialize_materials(const Mesh& mesh, const Material& materials)
{
    Materials mats;
    int num_cells = mesh.num_cells;
    mats.sigma_a.resize(num_cells);
    mats.heat_capacity.resize(num_cells);

    for (int i = 0; i < num_cells; ++i) {
        int material_id = mesh.cells[i].material_id;
        mats.sigma_a[i] = materials.get_sigma_a(material_id);
        mats.heat_capacity[i] = materials.get_heat_capacity(material_id);
    }

    return mats;
}*/


double absorption_coeff(int i, int j) {
    double siga;
    double sigaupper=2000.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        siga = sigaupper ;
    } else {
        siga = 20.0;  //20.0
    }
    
    return siga;
}

double specific_heat(int i, int j) {
    double ca;
    double caupper=10000.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        ca = caupper ;
    } else {
        ca = 1;
    }
    
    return ca;
}

double initial_temperature(int i, int j) {
    double tini=300;
    if(i<4 &&  j< 3*(NY/2)/4   && j>NY/4) {
    //if(i<4 &&  j< (NY/2)/2   ) {
        tini=TINI; // top hot configuration
    } else {
        tini=300.0;
    }

    return tini;
}

Materials::Materials() {
        // Initialize material properties for different material IDs
        // For simplicity, we define properties for two materials: 0 and 1
        material_properties.clear();
        // Define properties for material ID 0 (e.g., void or air)
        material_properties[0] = {0.0, 1.0}; // sigma_a, heat_capacity

        // Define properties for material ID 1 (e.g., pipe material)
        material_properties[1] = {100.0, 10.0}; // sigma_a, heat_capacity
    }


    double Materials::get_sigma_a(int material_id)  {
        auto it = material_properties.find(material_id);
        if (it != material_properties.end()) {
            return it->second.sigma_a;
        }
        return 0.0; // Default value if material ID not found
    }

    double Materials::get_heat_capacity(int material_id)  {
        auto it = material_properties.find(material_id);
        if (it != material_properties.end()) {
            return it->second.heat_capacity;
        }
        return 1.0; // Default value if material ID not found
    }

    int Materials::set_sigma_a(int material_id, double siga)  {
        int status=0;
        auto it = material_properties.find(material_id);
        if (it != material_properties.end() && siga >= 0.0) {
            it->second.sigma_a=siga;
            status=1.0;
        }
        return status; // Default value if material ID not found
    }

    int Materials::set_heat_capacity(int material_id, double heatcapacity)  {
        int status=0;
        auto it = material_properties.find(material_id);
        if (it != material_properties.end()) {
            it->second.heat_capacity=heatcapacity;
            status=1.0;
        }
        return status; // Default value if material ID not found
    }   