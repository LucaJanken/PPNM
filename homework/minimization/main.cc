#include <iostream>
#include "min.h"

// Define Rosenbrock's valley function
double rosenbrock(const vector& x) {
    return std::pow(1 - x[0], 2) + 100 * std::pow(x[1] - std::pow(x[0], 2), 2);
}

// Define Himmelblau's function
double himmelblau(const vector& x) {
    return std::pow(x[0] * x[0] + x[1] - 11, 2) + std::pow(x[0] + x[1] * x[1] - 7, 2);
}

// Define Breit-Wigner deviation function
double deviation(const vector& params, const vector& E, const vector& sigma, const vector& delta_sigma) {
    double m = params[0];
    double Gamma = params[1];
    double A = params[2];
    double sum = 0;
    for (size_t i = 0; i < E.size; ++i) {
        double term = A / ((E[i] - m) * (E[i] - m) + (Gamma * Gamma) / 4);
        sum += std::pow((term - sigma[i]) / delta_sigma[i], 2);
    }
    return sum;
}

int main() {
    // Set accuracy goal
    double tol = 1e-3;

    // Rosenbrock's valley function
    vector ros_init(2);
    ros_init[0] = 0;
    ros_init[1] = 0;
    std::cout << "Starting Rosenbrock's valley optimization..." << std::endl;
    std::cout << "Expect minimum at (1, 1)" << std::endl;
    vector min_rosenbrock = quasi_newton(rosenbrock, ros_init, tol, false);
    std::cout << "Found minimum of Rosenbrock's valley function with (x0, y0) = (0.0, 0.0) at: "
              << "(" << min_rosenbrock[0] << ", " << min_rosenbrock[1] << ")" << std::endl;
    
    ros_init[0] = 2.0;
    ros_init[1] = 0;
    vector min_rosenbrock2 = quasi_newton(rosenbrock, ros_init, tol, false);
    std::cout << "Found minimum of Rosenbrock's valley function with (x0, y0) = (2.0, 0.0) at: "
              << "(" << min_rosenbrock2[0] << ", " << min_rosenbrock2[1] << ")" << std::endl;

    ros_init[0] = 0;
    ros_init[1] = 2.0;
    vector min_rosenbrock3 = quasi_newton(rosenbrock, ros_init, tol, false);
    std::cout << "Found minimum of Rosenbrock's valley function with (x0, y0) = (0.0, 2.0) at: "
              << "(" << min_rosenbrock3[0] << ", " << min_rosenbrock3[1] << ")" << std::endl;
    
    ros_init[0] = 2.0;
    ros_init[1] = 2.0;
    vector min_rosenbrock4 = quasi_newton(rosenbrock, ros_init, tol, false);
    std::cout << "Found minimum of Rosenbrock's valley function with (x0, y0) = (2.0, 2.0) at: "
              << "(" << min_rosenbrock4[0] << ", " << min_rosenbrock4[1] << ")" << std::endl;

    // Himmelblau's function
    vector him_init(2);
    him_init[0] = -2.0;
    him_init[1] = -2;
    std::cout << "Starting Himmelblau's function optimization..." << std::endl;
    std::cout << "Expect minima at (-3.779310, -3.283186), (3.584428, -1.848126), (3.0, 2.0), (-2.805118, 3.131312)" << std::endl;
    vector min_himmelblau = quasi_newton(himmelblau, him_init, tol, false);
    std::cout << "Found minimum of Himmelblau's function with (x0, y0) = (-2.0, -2.0) at: "
              << "(" << min_himmelblau[0] << ", " << min_himmelblau[1] << ")" << std::endl;

    him_init[0] = 4.0;
    him_init[1] = -2.0;
    vector min_himmelblau2 = quasi_newton(himmelblau, him_init, tol, false);
    std::cout << "Found minimum of Himmelblau's function with (x0, y0) = (4.0, -2.0) at: "
              << "(" << min_himmelblau2[0] << ", " << min_himmelblau2[1] << ")" << std::endl;
    
    him_init[0] = 4.0;
    him_init[1] = 2.0;
    vector min_himmelblau3 = quasi_newton(himmelblau, him_init, tol, false);
    std::cout << "Found minimum of Himmelblau's function with (x0, y0) = (4.0, 2.0) at: "
              << "(" << min_himmelblau3[0] << ", " << min_himmelblau3[1] << ")" << std::endl;

    him_init[0] = -3.0;
    him_init[1] = 4.0;
    vector min_himmelblau4 = quasi_newton(himmelblau, him_init, tol, false);
    std::cout << "Found minimum of Himmelblau's function with (x0, y0) = (-3.0, 4.0) at: "
              << "(" << min_himmelblau4[0] << ", " << min_himmelblau4[1] << ")" << std::endl;

    // Task 2
    vector energy, signal, error;
    double e, s, err;

    while (std::cin >> e >> s >> err) {
        energy.push_back(e);
        signal.push_back(s);
        error.push_back(err);
    }

    vector breit_init(3);
    breit_init[0] = 125;
    breit_init[1] = 2;
    breit_init[2] = 8;
    
    vector params = quasi_newton([&](const vector& v) { return deviation(v, energy, signal, error); }, breit_init, tol, false);
    
    std::cout << "Fitting parameters (m, Γ, A):" << std::endl;
    std::cout << params[0] << " " << params[1] << " " << params[2] << std::endl;

    // Task 3
    vector params2 = quasi_newton([&](const vector& v) { return deviation(v, energy, signal, error); }, breit_init, tol, true);

    std::cout << "Fitting parameters (m, Γ, A) with SR1 update:" << std::endl;
    std::cout << params2[0] << " " << params2[1] << " " << params2[2] << std::endl;
    
    return 0;
}
