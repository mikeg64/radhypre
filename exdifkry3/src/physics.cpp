#include "../include/physics.h"

 

void linearize_emissive_source(PhysicsState& state) {
    for (int i = 0; i < state.mesh.num_cells; ++i) {
        double T = state.temperature[i];
        double T4 = std::pow(T, 4);
        for (int g = 0; g < NUM_GROUPS; ++g) {
            state.source_term[g][i] = state.sigma_a[g][i] * STEFAN_BOLTZMANN * T4;
        }
    }
}

 

void solve_material_heating(const Mesh& mesh, PhysicsState& state) {
    for (int i = 0; i < mesh.num_cells; ++i) {
        double delta_E = 0.0;
        for (int g = 0; g < NUM_GROUPS; ++g) {
            delta_E += state.sigma_a[g][i] * (state.radiation_flux[g][i] - state.source_term[g][i]);
        }
        state.temperature[i] += DT * delta_E / state.heat_capacity[i];
    }

}

State initialize_physics(const Mesh& mesh, const Materials& materials) {

    State state;
    int num_cells = mesh.num_cells;

    // Resize group-dependent vectors
    state.radiation_flux.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
    state.sigma_a.resize(NUM_GROUPS, std::vector(num_cells, 0.0));
    state.source_term.resize(NUM_GROUPS, std::vector(num_cells, 0.0));

    // Resize temperature and heat capacity
    state.temperature.resize(num_cells, 0.0);
    state.heat_capacity.resize(num_cells, 0.0);

    for (int i = 0; i < num_cells; ++i) {
        int mat_id = mesh.cells[i].material_id;
        // Set initial temperature (e.g., 1 keV blackbody)
        state.temperature[i] = TINI; //1.16e7; // Kelvin
        // Set heat capacity from material database
        state.heat_capacity[i] = materials.heat_capacity[mat_id];
        // Set absorption coefficients per group
        for (int g = 0; g < NUM_GROUPS; ++g) {
            state.sigma_a[g][i] = materials.sigma_a[mat_id][g];
        }

        // Initialize source term using blackbody emission
        double T4 = std::pow(state.temperature[i], 4);
        for (int g = 0; g < NUM_GROUPS; ++g) {
            state.source_term[g][i] = state.sigma_a[g][i] * STEFAN_BOLTZMANN * T4;
        }
    }

    return state;

}

