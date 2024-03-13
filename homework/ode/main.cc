#include "rkode.h"
#include <iostream>
#include <fstream>

// Vector that represents the system of ODEs (Harmonic oscillator)
vector f(double x, const vector& y) {
    vector dy(2);
    dy[0] = y[1];        // u' = v
    dy[1] = -y[0];       // v' = -u
    return dy;
}

// Vector that represents the system of ODEs (Oscillator with friction)
vector pendulum(double t, const vector& y) {
    double b = 0.25; // Damping coefficient
    double c = 5.0; // Gravitational constant * length of pendulum (or similar constant)
    vector dy(2);
    dy[0] = y[1]; // theta'(t) = omega(t)
    dy[1] = -b*y[1] - c*sin(y[0]); // omega'(t) = -b*omega(t) - c*sin(theta(t))
    return dy;
}

// Vector that represents the system of ODEs (Equitorial motion)
vector equatorial_motion(double phi, const vector& y, double epsilon) {
    vector dy(2);
    dy[0] = y[1]; // y0' = y1
    dy[1] = -y[0] + 1 + epsilon * y[0] * y[0]; // y1' = -y0 + 1 + epsilon * y0^2
    return dy;
}

// Helper function to compute acceleration of one body due to the other two
vector acceleration(const vector& r, const vector& r1, const vector& r2) {
    vector diff1 = r1 - r;
    vector diff2 = r2 - r;
    return diff1 / pow(diff1.norm(), 3) + diff2 / pow(diff2.norm(), 3);
}

// Vector that represents the system of ODEs (figure-8 orbit)
vector figure8_equations(double t, const vector& z) {
    vector r1(2), r2(2), r3(2), v1(2), v2(2), v3(2);
    r1[0] = z[6]; r1[1] = z[7];
    r2[0] = z[8]; r2[1] = z[9];
    r3[0] = z[10]; r3[1] = z[11];
    v1[0] = z[0]; v1[1] = z[1];
    v2[0] = z[2]; v2[1] = z[3];
    v3[0] = z[4]; v3[1] = z[5];

    // Compute accelerations using Newton's law of universal gravitation
    vector a1 = acceleration(r1, r2, r3);
    vector a2 = acceleration(r2, r1, r3);
    vector a3 = acceleration(r3, r1, r2);

    // Construct z'
    vector z_prime(12);
    z_prime[0] = a1[0]; z_prime[1] = a1[1];
    z_prime[2] = a2[0]; z_prime[3] = a2[1];
    z_prime[4] = a3[0]; z_prime[5] = a3[1];
    z_prime[6] = v1[0]; z_prime[7] = v1[1];
    z_prime[8] = v2[0]; z_prime[9] = v2[1];
    z_prime[10] = v3[0]; z_prime[11] = v3[1];

    return z_prime;
}


int main() {
    
    // Example u'' = -u
    vector ystart(2);
    ystart[0] = 1;       // u(0) = 1
    ystart[1] = 0;       // u'(0) = 0

    std::pair<vector, std::vector<vector>> result = driver(f, {0, 2 * M_PI}, ystart, 0.1, 0.01, 0.01);

    std::ofstream shm_file("shm.csv");
    shm_file << "x u v\n";
    for (size_t i = 0; i < result.first.size; ++i) {
        shm_file << result.first[i] << " " << result.second[i][0] << " " << result.second[i][1] << "\n";
    }
    shm_file.close();

    // Example theta''(t) = -b*theta'(t) - c*sin(theta(t))
    vector ystart2(2);
    ystart2[0] = M_PI - 0.1; // theta(0) = pi - 0.1, nearly vertical
    ystart2[1] = 0; // omega(0) = 0, initially at rest

    std::pair<vector, std::vector<vector>> result2 = driver(pendulum, {0, 10}, ystart2, 0.1, 0.01, 0.01);

    std::ofstream shm_friction_file("shm_friction.csv");
    shm_friction_file << "t theta omega\n";
    for (size_t i = 0; i < result2.first.size; ++i) {
        shm_friction_file << result2.first[i] << " " << result2.second[i][0] << " " << result2.second[i][1] << "\n";
    }
    shm_friction_file.close();
    
    // Example of Newtonian Circular Motion (ε=0)
    double epsilonA = 0;
    vector ystartA(2);
    ystartA[0] = 1.0;  // u(0) = 1
    ystartA[1] = 0.0;  // u'(0) = 0
    std::function<vector(double)> driver_interpA = driver_interp([&](double phi, const vector& y){ return equatorial_motion(phi, y, epsilonA); },
                              {0, 2 * M_PI}, // Integration interval
                              ystartA); // Parameters: h, acc, eps
    
    std::ofstream ncm_file("ncm.csv");
    for (double x = 0; x <= 2 * M_PI; x += M_PI / 32) {
        auto y_interp = driver_interpA(x);
        ncm_file << x << " " << y_interp[0] << "\n";
    }

    ncm_file.close();


    // Example of Newtonian Elliptical Motion (ε=0)
    double epsilonB = 0;
    vector ystartB(2);
    ystartB[0] = 1.0;  // u(0) = 1
    ystartB[1] = -0.5; // u'(0) ≈ -0.5
    std::function<vector(double)> driver_interpB = driver_interp([&](double phi, const vector& y){ return equatorial_motion(phi, y, epsilonB); },
                         {0, 2 * M_PI}, // Extend the integration interval if necessary
                         ystartB);
    
    std::ofstream nem_file("nem.csv");
    double x = 0; // Declare the variable "x" before the try block
    try {
        for (x = 0; x <= 2 * M_PI; x += M_PI / 64) {
            auto y_interp = driver_interpB(x);
            nem_file << x << " " << y_interp[0] << "\n";
        }
    } catch (const std::invalid_argument& e) {
        std::cerr << "Interpolation error at phi = " << x << ": " << e.what() << std::endl;
    }

    nem_file.close();

    // Example Relativistic Precession (ε≈0.01)
    double epsilonC = 0.01;
    vector ystartC(2);
    ystartC[0] = 1.0;  // u(0) = 1
    ystartC[1] = -0.5; // u'(0) ≈ -0.5
    std::function<vector(double)> driver_interpC = driver_interp([&](double phi, const vector& y){ return equatorial_motion(phi, y, epsilonC); },
                         {0, 6 * M_PI}, // Extend the integration interval if necessary
                         ystartC);

    std::ofstream rp_file("rp.csv");
    for (double x = 0; x <= 6 * M_PI; x += M_PI / 32) {
        auto y_interp = driver_interpC(x);
        rp_file << x << " " << y_interp[0] << "\n";
    }

    rp_file.close();

    // Example of the figure 8 orbit
    vector z(12);
    
    z[0] = 0.466203685; z[1] = 0.43236573;   // v1
    z[2] = 0.466203685; z[3] = 0.43236573;   // v2
    z[4] = -0.93240737; z[5] = -0.86473146;  // v3
    z[6] = 0.97000436;  z[7] = -0.24308753;  // r1
    z[8] = -0.97000436; z[9] = 0.24308753;   // r2
    z[10] = 0;          z[11] = 0;           // r3

    std::function<vector(double)> figure8_interp = driver_interp(
        figure8_equations, {0, 6.32591398 / 3}, z, 0.01, 0.01, 0.01, 10000
    );

    std::ofstream fig8_file("figure8.csv");
    for (double t = 0; t <= 6.32591398 / 3; t += 0.01) {
        auto z_interp = figure8_interp(t);
        fig8_file << z_interp[6] << " " << z_interp[7] << ", "  // Position of body 1
                << z_interp[8] << " " << z_interp[9] << ", "  // Position of body 2
                << z_interp[10] << " " << z_interp[11] << "\n";  // Position of body 3
    }
    fig8_file.close();

    
    return 0;
}