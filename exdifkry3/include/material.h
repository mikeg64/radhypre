#pragma once
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
    Materials() {
        // Initialize material properties for different material IDs
        // For simplicity, we define properties for two materials: 0 and 1
        material_properties.clear();
        // Define properties for material ID 0 (e.g., void or air)
        material_properties[0] = {0.0, 1.0}; // sigma_a, heat_capacity

        // Define properties for material ID 1 (e.g., pipe material)
        material_properties[1] = {100.0, 10.0}; // sigma_a, heat_capacity
    }

    double get_sigma_a(int material_id) const {
        auto it = material_properties.find(material_id);
        if (it != material_properties.end()) {
            return it->second.sigma_a;
        }
        return 0.0; // Default value if material ID not found
    }

    double get_heat_capacity(int material_id) const {
        auto it = material_properties.find(material_id);
        if (it != material_properties.end()) {
            return it->second.heat_capacity;
        }
        return 1.0; // Default value if material ID not found
    }
};

double absorption_coeff(int i, int j) {
    double siga;
    double sigaupper=2000.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        siga = sigaupper ;
    } else {
        siga = 20.0;  //20.0
    }
    
    return siga;
};

double specific_heat(int i, int j) {
    double ca;
    double caupper=10000.0;

    if ((i < NX / 5   && j > NY / 2)   ||  (i > (4*NX / 5)   && j > NY / 2)   ||  (i > (2*NX / 5) &&  i<(3*NX/5)  && j < NY / 2) ) {
        ca = caupper ;
    } else {
        ca = 1;
    }
    
    return ca;
};

double initial_temperature(int i, int j) {
    double tini=300;
    if(i<4 &&  j< 3*(NY/2)/4   && j>NY/4) {
    //if(i<4 &&  j< (NY/2)/2   ) {
        tini=TINI; // top hot configuration
    } else {
        tini=300.0;
    }

    return tini;
};

Materials initialize_materials(const Mesh& mesh, const Material& materials);
Materials initialize_materials(const Mesh& mesh);
