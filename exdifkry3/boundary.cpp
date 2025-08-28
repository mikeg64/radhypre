#include "include/boundary.h"

#include <cmath>

 

// Milne boundary condition: compute reflected flux

double compute_milne_reflected_flux(double incident_flux, double temperature, double sigma_a) {

    // Milne extrapolation factor (approximate)

    const double milne_factor = 0.7104;

    double blackbody_flux = sigma_a * STEFAN_BOLTZMANN * std::pow(temperature, 4);

    return milne_factor * blackbody_flux + (1 - milne_factor) * incident_flux;

}

 

void apply_milne_boundary_conditions(const Mesh& mesh, PhysicsState& state) {

    for (auto& boundary : mesh.boundaries) {

        if (boundary.type == WALL) {

            for (int g = 0; g < NUM_GROUPS; ++g) {

                double incident_flux = state.radiation_flux[g][boundary.cell];

                double reflected_flux = compute_milne_reflected_flux(

                    incident_flux,

                    state.temperature[boundary.cell],

                    state.sigma_a[g][boundary.cell]

                );

                state.radiation_flux[g][boundary.cell] = reflected_flux;

            }

        }

    }

}