#include "roots.h"
#include "../ode/rkode.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

// Test 1: One-dimensional function, f(x) = x^2 - 2
vector testFunc1D(vector x) {
    vector f(1);
    f[0] = x[0]*x[0] - 2;
    return f;
}

// Test 2: Simple two-dimensional linear equations
vector testFunc2D(vector x) {
    vector f(2);
    f[0] = 2*x[0] + 3*x[1] - 5; // 2x1 + 3x2 = 5
    f[1] = 3*x[0] - x[1] - 2;   // 3x1 - x2 = 2
    return f;
}

// Test 3: Gradient of the Rosenbrock function
vector rosenbrockGrad(vector x) {
    vector grad(2);
    grad[0] = -2*(1 - x[0]) - 400*x[0]*(x[1] - x[0]*x[0]);
    grad[1] = 200*(x[1] - x[0]*x[0]);
    return grad;
}

// System of ODEs representing the Schrödinger equation for the hydrogen atom
vector schrodinger(double r, const vector& y, double E) {
    vector dy(2);
    dy[0] = y[1]; // f' = g
    dy[1] = -2.0 * (E + 1.0 / r) * y[0]; // g' = -2(E + 1/r)f
    return dy;
}

// Function to find f(rmax) for a given energy E
double M(double E, double rmin, double rmax, const vector& ystart, double acc, double eps) {
    auto result = driver([&](double r, const vector& y){ return schrodinger(r, y, E); }, 
                         {rmin, rmax}, ystart, acc, eps);
    return result.second.back()[0];
}

// Wrapper for M(E) to use with Newton's method
vector ME_wrapper(vector E, double rmin, double rmax, double acc, double eps) {
    vector ystart(2);
    ystart[0] = rmin - pow(rmin, 2); 
    ystart[1] = 1 - 2 * rmin;    
    

    double valueAtRmax = M(E[0], rmin, rmax, ystart, acc, eps);
    vector result(1);
    result[0] = valueAtRmax; 
    return result;
}

double E0Param(double rmin, double rmax, double acc, double eps) {
    vector initialGuess(1);
    initialGuess[0] = -1; 
    
    vector E0 = newton([&](vector E){ return ME_wrapper(E, rmin, rmax, acc, eps); }, initialGuess, eps);
    return E0[0]; 
}

// Function to match the numerical solution to the asymptotic form
double matchAsymptoticForm(double E, double r, const vector& y) {
    double k = std::sqrt(-2 * E);
    double f_asymptotic = r * std::exp(-k * r);
    double df_asymptotic = (1 - k * r) * std::exp(-k * r);
    double delta_f = y[0] - f_asymptotic;
    double delta_df = y[1] - df_asymptotic;

    return std::abs(delta_f) + std::abs(delta_df);
}

// Modified function to find f(rmax) that incorporates the improved boundary condition
double M_mod(double E, double rmin, double rmax, const vector& ystart, double acc, double eps) {
    auto result = driver([&](double r, const vector& y) { return schrodinger(r, y, E); },
                         {rmin, rmax}, ystart, acc, eps);
    double discrepancy = matchAsymptoticForm(E, rmax, result.second.back());
    return discrepancy;
}

double ImprovedE0(double rmin, double rmax, double acc, double eps, double EnergyGuess = -1.0) {
    vector initialGuess(1);
    initialGuess[0] = EnergyGuess;

    auto energyFinder = [&](vector E) -> vector {
        double k = std::sqrt(-2 * E[0]);
        vector ystart(2);
        ystart[0] = rmin * std::exp(-k * rmin);
        ystart[1] = (1 - k * rmin) * std::exp(-k * rmin); 

        auto result = driver([&](double r, const vector& y) { return schrodinger(r, y, E[0]); },
                             {rmin, rmax}, ystart, acc, eps);
        double discrepancy = matchAsymptoticForm(E[0], rmax, result.second.back());

        vector discrepancyVec(1);
        discrepancyVec[0] = discrepancy; 
        return discrepancyVec;
    };
    
    vector foundE = newton(energyFinder, initialGuess, eps);
    return foundE[0]; 
}

int main() {
    vector initialGuess, root;

    // Test 1: One-dimensional function
    initialGuess = vector(1);
    std::cout << "Testing simple one-dimensional function x^2 - 2" << std::endl;
    initialGuess[0] = -2; // Starting guess
    root = newton(testFunc1D, initialGuess);
    std::cout << "Initial guess: " << initialGuess[0] << ", Solution found at: " << root[0] << std::endl;
    initialGuess[0] = -1; // Starting guess
    root = newton(testFunc1D, initialGuess);
    std::cout << "Initial guess: " << initialGuess[0] << ", Solution found at: " << root[0] << std::endl;
    initialGuess[0] = 1; // Starting guess
    root = newton(testFunc1D, initialGuess);
    std::cout << "Initial guess: " << initialGuess[0] << ", Solution found at: " << root[0] << std::endl;
    initialGuess[0] = 2; // Starting guess
    root = newton(testFunc1D, initialGuess);
    std::cout << "Initial guess: " << initialGuess[0] << ", Solution found at: " << root[0] << std::endl;
    
    // Test 2: Simple two-dimensional linear equations
    std::cout << "Testing simple two-dimensional function 2x1 + 3x2 = 5, 3x1 - x2 = 2" << std::endl;
    initialGuess = vector(2);
    initialGuess[0] = 0; // Starting guess for x1
    initialGuess[1] = 0; // Starting guess for x2
    root = newton(testFunc2D, initialGuess);
    std::cout << "Initial guess: (" << initialGuess[0] << ", " << initialGuess[1] << ")" << ", Solution found at: (" << root[0] << ", " << root[1] << ")" << std::endl;
    initialGuess[0] = 2; // Starting guess for x1
    initialGuess[1] = 0; // Starting guess for x2
    root = newton(testFunc2D, initialGuess);
    std::cout << "Initial guess: (" << initialGuess[0] << ", " << initialGuess[1] << ")" << ", Solution found at: (" << root[0] << ", " << root[1] << ")" << std::endl;
    initialGuess[0] = 0; // Starting guess for x1
    initialGuess[1] = 2; // Starting guess for x2
    root = newton(testFunc2D, initialGuess);
    std::cout << "Initial guess: (" << initialGuess[0] << ", " << initialGuess[1] << ")" << ", Solution found at: (" << root[0] << ", " << root[1] << ")" << std::endl;
    initialGuess[0] = 2; // Starting guess for x1
    initialGuess[1] = 2; // Starting guess for x2
    root = newton(testFunc2D, initialGuess);
    std::cout << "Initial guess: (" << initialGuess[0] << ", " << initialGuess[1] << ")" << ", Solution found at: (" << root[0] << ", " << root[1] << ")" << std::endl;
    
    // Test 3: Rosenbrock's valley function
    std::cout << "Testing Rosenbrock function..." << std::endl;
    initialGuess = vector(2);
    initialGuess[0] = 0; // Starting guess for x1
    initialGuess[1] = 0;  // Starting guess for x2
    root = newton(rosenbrockGrad, initialGuess);
    std::cout << "Initial guess: (" << initialGuess[0] << ", " << initialGuess[1] << ")" << ", Solution found at: (" << root[0] << ", " << root[1] << ")" << std::endl;
    initialGuess[0] = 2; // Starting guess for x1
    initialGuess[1] = 0;  // Starting guess for x2
    root = newton(rosenbrockGrad, initialGuess);
    std::cout << "Initial guess: (" << initialGuess[0] << ", " << initialGuess[1] << ")" << ", Solution found at: (" << root[0] << ", " << root[1] << ")" << std::endl;
    initialGuess[0] = 0; // Starting guess for x1
    initialGuess[1] = 2;  // Starting guess for x2
    root = newton(rosenbrockGrad, initialGuess);
    std::cout << "Initial guess: (" << initialGuess[0] << ", " << initialGuess[1] << ")" << ", Solution found at: (" << root[0] << ", " << root[1] << ")" << std::endl;
    initialGuess[0] = 2; // Starting guess for x1
    initialGuess[1] = 2;  // Starting guess for x2
    root = newton(rosenbrockGrad, initialGuess);
    std::cout << "Initial guess: (" << initialGuess[0] << ", " << initialGuess[1] << ")" << ", Solution found at: (" << root[0] << ", " << root[1] << ")" << std::endl;

    // Part B

    // Initial guess for E
    initialGuess = vector(1);
    initialGuess[0] = -1;

    // Using Newton's method to find E0, the root of M(E) = 0
    vector E0 = newton([&](vector E){ return ME_wrapper(E, 1e-4, 8.0, 0.01, 0.01); }, initialGuess);
    std::cout << "Lowest energy level found: E0 = " << E0[0] << std::endl;

    // Solve the ODE with the found E0 to obtain the wave function
    double rmin = 1e-4, rmax = 8.0;
    vector ystart(2);
    ystart[0] = rmin - pow(rmin, 2); // f(rmin) = rmin - rmin^2
    ystart[1] = 1 - 2*rmin; // This is a derivative approximation of f at rmin
    auto result = driver([&](double r, const vector& y) { return schrodinger(r, y, E0[0]); },
                         {rmin, rmax}, ystart);

    // Write the wave function to a CSV file for plotting
    std::ofstream wavefunction_file("wavefunction.csv");
    for (size_t i = 0; i < result.first.size; ++i) {
        wavefunction_file << result.first[i] << " " << result.second[i][0] << "\n";
    }
    wavefunction_file.close();

    // Convergence test
    int N = 10;
    vector rmins(N), rmaxs(N), accs(N), epss(N);
    vector Ermin(N), Ermax(N), Eacc(N), Eeps(N);

    for (int i = 0; i < N; ++i) {
        rmins[i] = 1 / std::pow(2, i);
        rmaxs[i] = 1 + i;
        accs[i] = 500 / std::pow(10, i);
        epss[i] = 500 / std::pow(10, i);
    }

    for (int i = 0; i < N; ++i) {
        Ermin[i] = E0Param(rmins[i], 8.0, 0.01, 0.01);
    }

    for (int i = 0; i < N; ++i) {
        Ermax[i] = E0Param(1e-4, rmaxs[i], 0.01, 0.01);
    }

    for (int i = 0; i < N; ++i) {
        Eacc[i] = E0Param(1e-4, 8.0, accs[i], 1);
    }

    for (int i = 0; i < N; ++i) {
        Eeps[i] = E0Param(1e-4, 8.0, 0.01, epss[i]);
    }

    std::ofstream out("convergence.csv");
    out << "rmin Ermin rmax Ermax acc Eacc eps Eeps\n";
    out << std::fixed << std::setprecision(10);
    for (int i = 0; i < N; ++i) {
        out << rmins[i] << " " << Ermin[i] << " "
            << rmaxs[i] << " " << Ermax[i] << " "
            << accs[i] << " " << Eacc[i] << " "
            << epss[i] << " " << Eeps[i] << "\n";
    }

    out.close();

    // Part C

    // Test convergence with improved boundary condition
    vector E0s(N);
    for (int i = 0; i < N; ++i) {
        E0s[i] = ImprovedE0(1e-4, 1 + i, 0.01, 0.01);
    }

    std::ofstream out2("convergence_improved.csv");
    out2 << "rmax E0\n";
    out2 << std::fixed << std::setprecision(10);
    for (int i = 0; i < N; ++i) {
        out2 << 1 + i << " " << E0s[i] << "\n";
    }

    out2.close();

    return 0;
}