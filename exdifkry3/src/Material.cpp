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
