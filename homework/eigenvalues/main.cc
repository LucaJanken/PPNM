#include "JDWCS.h"
#include <random>
#include <iostream>
#include <fstream>
#include <cstring>
#include <chrono>

// Function to write data to CSV
void write_to_csv(const std::string &filename, const vector &x_values, const vector &y_values) {
    std::ofstream file(filename);
    for (size_t i = 0; i < x_values.size; ++i) {
        file << x_values[i] << " " << y_values[i] << std::endl;
    }
    file.close();
}

// Function to save the first few eigenfunctions to CSV
void save_eigenfunctions(const matrix& eig_vecs, double dr, int rmax, int number = 5) {
    std::ofstream file("eigenfunctions.csv"); // Fixed filename
    int npoints = static_cast<int>(rmax / dr) - 1;
    double Const = 1.0 / sqrt(dr); // Normalization constant

    for (int i = 0; i < npoints; ++i) {
        double r = dr * (i + 1); // Radial position
        file << r; // Write the radial position first

        for (int k = 0; k < number; ++k) {
            double normalized_value = Const * eig_vecs.get(i, k); // Normalize the eigenfunction value
            file << " " << normalized_value; // Append eigenfunction values separated by spaces
        }
        file << std::endl; // Newline after each row
    }

    file.close();
}

// Function to generate a random symmetric matrix
matrix rnd_symmetric_matrix(size_t n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    matrix A(n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n; j++) {
            double value = dis(gen);
            A.set(i, j, value);
            if (i != j) {
                A.set(j, i, value);
            }
        }
    }
    return A;
}

// Function to generate Hamiltonian matrix
matrix Hamiltonian(int rmax, double dr) {
    int npoints = static_cast<int>(rmax / dr) - 1;
    vector r(npoints);
    for (int i = 0; i < npoints; i++) {
        r[i] = dr * (i + 1);
    }
    matrix H(npoints, npoints);
    for (int i = 0; i < npoints - 1; i++) {
        H.set(i, i, -2);
        H.set(i, i + 1, 1);
        H.set(i + 1, i, 1); 
    }
    H.set(npoints - 1, npoints - 1, -2);
    H = H * (-0.5 / dr / dr);
    for (int i = 0; i < npoints; i++) {
        H.set(i, i, H.get(i, i) - 1 / r[i]);
    }
    return H;
}

int main(int argc, char* argv[]) {
    double dr = 0.05; // Default value
    int rmax = 25; // Default value

    // Parse command-line arguments
    for(int i = 1; i < argc; ++i) {
        if(std::strcmp(argv[i], "-dr") == 0) {
            if(i + 1 < argc) { // Make sure we don't go out of bounds
                dr = std::atof(argv[++i]); // Increment i and convert the next argument to double
            }
        } else if(std::strcmp(argv[i], "-rmax") == 0) {
            if(i + 1 < argc) { // Make sure we don't go out of bounds
                rmax = std::atoi(argv[++i]); // Increment i and convert the next argument to int
            }
        }
    }

    // Generate a random symmetric matrix
    int m = 5;
    matrix A = rnd_symmetric_matrix(m);
    A.print("A:");

    // Perform JDWCS
    std::pair<matrix, matrix> result = JDWCS::cyclic(A);
    matrix D = result.first;
    matrix V = result.second;

    D.print("D:");
    V.print("V:");

    // Verify V^T * A * V = D
    matrix VT = V.transpose();
    matrix VTA = VT * A; 
    matrix VTAV = VTA * V;

    VTAV.print("V^T * A * V:");
    if (matrix::compare(VTAV, D)) {
        std::cout << "V^T * A * V is approximately equal to D.\n";
    } else {
        std::cout << "V^T * A * V is not approximately equal to D.\n";
    }

    // Verify V * D * V^T = A
    matrix VD = V * D; 
    matrix VDVT = VD * VT;

    VDVT.print("V * D * V^T:");
    if (matrix::compare(VDVT, A)) {
        std::cout << "V * D * V^T is approximately equal to A.\n";
    } else {
        std::cout << "V * D * V^T is not approximately equal to A.\n";
    }

    // Verify VT * V = I
    matrix VTV = VT * V; 
    VTV.print("V^T * V:");
    if (matrix::compare(VTV, matrix::identity(m))) {
        std::cout << "V^T * V is approximately equal to I.\n";
    } else {
        std::cout << "V^T * V is not approximately equal to I.\n";
    }

    // Verify V * VT = I
    matrix VVT = V * VT;
    VVT.print("V * V^T:");
    if (matrix::compare(VVT, matrix::identity(m))) {
        std::cout << "V * V^T is approximately equal to I.\n";
    } else {
        std::cout << "V * V^T is not approximately equal to I.\n";
    }

    // Generate Hamiltonian matrix
    matrix H = Hamiltonian(rmax, dr);

    // Perform JDWCS
    std::pair<matrix, matrix> result2 = JDWCS::cyclic(H);
    matrix eig_vals = result2.first;
    matrix eig_vecs = result2.second;

    // Save the first few eigenfunctions to CSV
    save_eigenfunctions(eig_vecs, dr, rmax);

    // Loop over dr and calculate ε0
    vector delta_r_values(10);
    for (size_t i = 0; i < 10; ++i) {
        double dr = 0.1 * (i + 1); 
        delta_r_values.set(i, dr); 
    }
    vector epsilon0_delta_r(10);

    for (size_t i = 0; i < delta_r_values.size; ++i) {
        matrix H_i = Hamiltonian(rmax, delta_r_values[i]);
        std::pair<matrix, matrix> result_i = JDWCS::cyclic(H_i);
        matrix eig_vals_i = result_i.first;
        epsilon0_delta_r.set(i, eig_vals_i.get(0, 0));
    }

    // Write ε0 vs Δr to CSV
    write_to_csv("delta_r_vs_epsilon0.csv", delta_r_values, epsilon0_delta_r);

    // Loop over rmax and calculate ε0
    vector rmax_values(10);
    for (size_t i = 0; i < 10; ++i) {
        int rmax = 1 + i; 
        rmax_values.set(i, rmax); 
    }
    vector epsilon0_rmax(10);

    for (size_t i = 0; i < rmax_values.size; ++i) {
        matrix H_i = Hamiltonian(rmax_values[i], dr);
        std::pair<matrix, matrix> result_i = JDWCS::cyclic(H_i);
        matrix eig_vals_i = result_i.first;
        epsilon0_rmax.set(i, eig_vals_i.get(0, 0));
    }

    // Write ε0 vs rmax to CSV
    write_to_csv("rmax_vs_epsilon0.csv", rmax_values, epsilon0_rmax);

    return 0;
}
