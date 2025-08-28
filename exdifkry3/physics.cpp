#include "include/physics.h"

 

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

