#pragma once
#include "setup.h"
#include "geometry.h"
#include <map>

class MaterialProperties {
    
public:
        double sigma_a;      // Absorption coefficient
        double heat_capacity; // Heat capacity
    };
// Setting up the material properties for the simulation
class Materials {
public:


    std::map<int, MaterialProperties> material_properties;
    Materials();

    double get_sigma_a(int material_id) ;

    double get_heat_capacity(int material_id);

    int set_sigma_a(int material_id, double siga);

    int set_heat_capacity(int material_id, double heatcapacity);

    int add_material(int material_id, const MaterialProperties& props);





};



//Materials initialize_materials(const Mesh& mesh, const Materials& materials);
Materials initialize_materials(const Mesh& mesh);
double absorption_coeff(int i, int j);
double specific_heat(int i, int j);
double initial_temperature(int i, int j);
