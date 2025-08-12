#pragma once



//same as B_nu but for Planck's law
// This function computes the spectral radiance of a black body at frequency nu and temperature T
double planck_emission(double nu, double T) {
    //const double h = 6.62607015e-27;     // erg·s
    //const double c = 2.99792458e10;      // cm/s
    //const double k = 1.380649e-16;       // erg/K

    double numerator = 2.0 * h * std::pow(nu, 3) / (c * c);
    double exponent = h * nu / (k * T);
    double denominator = std::exp(exponent) - 1.0;

    return numerator / denominator; // erg·s⁻¹·cm⁻²·Hz⁻¹·sr⁻¹
}


// Function to compute Planck distribution
double B_nu(double T, double nu, double dbnu=-1.0) {
    double tot=0;
    double nud=nu;

    if(dbnu==-1.0) {
        tot = (2 * h * std::pow(nu, 3) / std::pow(c, 2)) / (std::exp(h * nu / (k * T)) - 1);
    } else {
        nud=nu-(dbnu/2.0);
        nu=(nud>0?nud:nu); // Avoid division by zero
        tot = (2 * h * std::pow(nu, 3) / std::pow(c, 2)) / (std::exp(h * nu / (k * T)) - 1) * dbnu;
    }


    tot= (2 * h * pow(nu, 3) / pow(c, 2)) / (exp(h * nu / (k * T)) - 1);

    return tot;
}

// Function to compute Planck distribution
double dBnudT(double T, double nu) {
    //return 0.0;
    return (exp(h*nu/(k*T))*2 * pow(h,2) * pow(nu, 4) / (k*pow(T,2)*pow(c, 2))) / pow(exp(h * nu / (k * T)) - 1,2);
}


double dBdT(double nu, double T) {
    const double h = 6.62607015e-27;     // erg·s
    const double c = 2.99792458e10;      // cm/s
    const double k = 1.380649e-16;       // erg/K

    double x = h * nu / (k * T);
    double exp_x = std::exp(x);
    double factor = (2.0 * h * h * nu * nu * nu * nu) / (c * c * k * T * T);
    return factor * exp_x / std::pow(exp_x - 1.0, 2);
}


// Define the opacity function
std::function<double(double)> make_kappa_nu(double kappa0, double nu0, double alpha) {
    return [=](double nu) {
        return kappa0 * std::pow(nu / nu0, -alpha);
    };
}
// Function to compute the Rosseland mean opacity for a group of frequencies
// kappa_nu is a function that returns the opacity at frequency nu

//example usage of function
//    double T = 6000.0;                  // Temperature in Kelvin
//    double nu_min = 1e13;              // Minimum frequency in Hz
//    double nu_max = 1e15;              // Maximum frequency in Hz
//    double kappa0 = 1.0;               // Reference opacity
//    double nu0 = 1e14;                 // Reference frequency
//    double alpha = 1.5;                // Power-law index

    // Create the opacity function
//    auto kappa_nu = make_kappa_nu(kappa0, nu0, alpha);

    // Compute Rosseland mean opacity for the group
//    double kappa_R = rosseland_mean_opacity_group(kappa_nu, T, nu_min, nu_max);








double rosseland_mean_opacity_group(
    std::function<double(double)> kappa_nu, // frequency-dependent opacity
    double T,
    double nu_min,
    double nu_max,
    int N = 1000 // number of integration points
) {
    double dnu = (nu_max - nu_min) / N;
    double numerator = 0.0;
    double denominator = 0.0;

    for (int i = 0; i < N; ++i) {
        double nu = nu_min + i * dnu;
        double dB_dT = dBdT(nu, T);
        numerator += (1.0 / kappa_nu(nu)) * dB_dT * dnu;
        denominator += dB_dT * dnu;
    }

    return denominator / numerator; // Rosseland mean opacity for the group
}

double rosseland_mean_opacity_group(
    double kappa, // frequency-dependent opacity
    double T,
    int num_freq_bins = 1000 // number of integration points
) {
    double dnu = 1.0e14;
    double numerator = 0.0;
    double denominator = 0.0;

 


    for (int i = 0; i < num_freq_bins; ++i) {
        double nu = (i+1)*1e14;
        double dB_dT = dBdT(nu, T);
        numerator += (1.0 / kappa) * dB_dT * dnu;
        denominator += dB_dT * dnu;
    }

    return denominator / numerator; // Rosseland mn opacity for the group
}






double planck_source(double T) {
    return 4.0 * 5.67e-8 * std::pow(T, 4);  // Simplified gray approx.
}
