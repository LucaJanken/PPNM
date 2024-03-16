#include <iostream>
#include <cmath>
#include <functional>
#include <fstream>
#include "quadrature.h"


// Define the functions to integrate
double f_sqrt(double x) {
    return std::sqrt(x);
}

double f_inv_sqrt(double x) {
    return 1 / std::sqrt(x);
}

double f_cirle(double x) {
    return 4 * std::sqrt(1 - x * x);
}

double f_ln_sqrt(double x) {
    return std::log(x) / std::sqrt(x);
}

// Define function to calculate absolute error
double abs_error(double x, double y) {
    return std::abs(x - y);
}

// Error function
// For 0 <= z <= 1
double erfIntegral(double z, int& evaluations) {
    auto result = integrate([](double x) { return std::exp(-x*x); }, 0, z, 1e-8, 1e-8);
    evaluations += std::get<2>(result); // Update evaluations
    return std::get<0>(result); // Return the integral result
}

// For z > 1
double erfIntegralGreater(double z, int& evaluations) {
    auto result = integrate([z](double t) { return std::exp(-std::pow(z + (1 - t) / t, 2)) / (t * t); }, 0, 1, 1e-8, 1e-8);
    evaluations += std::get<2>(result); // Update evaluations
    return std::get<0>(result); // Return the integral result
}

// Definition of the error function
double erfFunction(double z, int& evaluations) {
    const double sqrtPi = std::sqrt(M_PI);
    evaluations = 0; // Reset or initialize evaluations count if necessary
    if (z < 0) {
        return -erfFunction(-z, evaluations);
    } else if (z <= 1) {
        return (2 / sqrtPi) * erfIntegral(z, evaluations);
    } else {
        return 1 - (2 / sqrtPi) * erfIntegralGreater(z, evaluations);
    }
}

// Approximation of the error function
double erf_approx(double x) {
    if (x < 0) return -erf(-x);
    double a[] = {0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429};
    double t = 1.0 / (1.0 + 0.3275911 * x);
    double sum = t * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4]))));
    return 1.0 - sum * exp(-x * x);
}

int main() {
    // Initialize evaluation counter
    int evaluations = 0;

    // Results
    auto result_sqrt = integrate(f_sqrt, 0, 1, 0.001, 0.001);
    auto result_inv_sqrt = integrate(f_inv_sqrt, 0, 1, 0.001, 0.001);
    auto result_circle = integrate(f_cirle, 0, 1, 0.001, 0.001);
    auto result_ln_sqrt = integrate(f_ln_sqrt, 0, 1, 0.001, 0.001);

    std::cout << "Integral of sqrt(x) over [0, 1]: " << std::get<0>(result_sqrt) << " Expected: 2/3" << " Error: " << abs_error(std::get<0>(result_sqrt), 2.0/3.0) << std::endl;
    std::cout << "Integral of 1/sqrt(x) over [0, 1]: " << std::get<0>(result_inv_sqrt) << " Expected: 2"  << " Error: " << abs_error(std::get<0>(result_inv_sqrt), 2.0) << std::endl;
    std::cout << "Integral of 4*sqrt(1-x^2) over [0, 1]: " << std::get<0>(result_circle) << " Expected: pi" << " Error: " << abs_error(std::get<0>(result_circle), M_PI) << std::endl;
    std::cout << "Integral of ln(x)/sqrt(x) over [0, 1]: " << std::get<0>(result_ln_sqrt) << " Expected: -4" << " Error: " << abs_error(std::get<0>(result_ln_sqrt), -4.0) << std::endl;

    // Calculate the error function for [-3, 3] and save data to erf_plot.csv file.
    std::ofstream file("erf_plot.csv");
    for (double z = -3; z <= 3; z += 0.1) {
        evaluations = 0; // Reset the counter for each calculation
        auto erfResult = erfFunction(z, evaluations); // Adjusted for tuple return
        file << z << " " << erfResult << std::endl;
    }
    file.close();

    // Compute differences, and write to CSV
    std::ifstream infile("../../exercises/plots/erf_tab.txt");
    std::ofstream outfile("diff.csv");

    double x, erf_tab, unused;

    while (infile >> x >> erf_tab >> unused) {
        evaluations = 0; // Reset the counter for each calculation
        auto erfResult = erfFunction(x, evaluations); // Assuming erfFunction now returns just the function value
        evaluations = 0; // Reset for approximation
        double diff_approx = std::abs(erf_approx(x) - erf_tab);
        outfile << x << " " << std::abs(erfResult - erf_tab) << " " << diff_approx << "\n";
    }

    infile.close();
    outfile.close();

    std::cout << "Plot and Comparison data of erf(x) written to erf_plot.csv and diff.csv" << std::endl;

    // Singular integrals and comparison
    evaluations = 0;
    auto result_inv_sqrt2 = integrate(f_inv_sqrt, 0, 1, 0.001, 0.001);
    auto result_inv_sqrt_trans = transint(f_inv_sqrt, 0, 1, 0.001, 0.001);
    std::cout << "Standard integration of 1/sqrt(x): Result = " << std::get<0>(result_inv_sqrt2) << " Expected: 2" << ", Evaluations = " << std::get<2>(result_inv_sqrt2) << std::endl;
    std::cout << "Clenshaw-Curtis transformation of 1/sqrt(x): Result = " << std::get<0>(result_inv_sqrt_trans) << " Expected: 2" << ", Evaluations = " << std::get<2>(result_inv_sqrt_trans) << std::endl;

    evaluations = 0;
    auto result_ln_sqrt3 = integrate(f_ln_sqrt, 0, 1, 0.001, 0.001);
    auto result_ln_sqrt_trans2 = transint(f_ln_sqrt, 0, 1, 0.001, 0.001);
    std::cout << "Standard integration of ln(x)/sqrt(x): Result = " << std::get<0>(result_ln_sqrt3) << " Expected: -4" << ", Evaluations = " << std::get<2>(result_ln_sqrt3) << std::endl;
    std::cout << "Clenshaw-Curtis transformation of ln(x)/sqrt(x): Result = " << std::get<0>(result_ln_sqrt_trans2) << " Expected: -4" << ", Evaluations = " << std::get<2>(result_ln_sqrt_trans2) << std::endl;

    // Test function: integral from 0 to infinity of exp(-x) dx = 1
    auto result_inf = integrate_infinite([](double x) { return std::exp(-x); }, 0, std::numeric_limits<double>::infinity(), 1e-8, 1e-8);
    std::cout << "exp(-x) from 0 to infinity: Result = " << std::get<0>(result_inf) << " Expected: 1" << ", Evaluations = " << std::get<2>(result_inf) << ", Error = " << abs_error(std::get<0>(result_inf), 1.0) << std::endl;

    // Test function: integral from -infinity to infinity of exp(-x^2) dx = sqrt(pi)
    auto result_inf2 = integrate_infinite([](double x) { return std::exp(-x*x); }, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 1e-8, 1e-8);
    std::cout << "exp(-x^2) from -infinity to infinity: Result = " << std::get<0>(result_inf2) << " Expected: sqrt(pi) = 1.772453851" << ", Evaluations = " << std::get<2>(result_inf2) << ", Error = " << abs_error(std::get<0>(result_inf2), std::sqrt(M_PI)) << std::endl;
    
    return 0;
}