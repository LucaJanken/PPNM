#include "ScalingInvestigation.h"
#include "EigenSolver.h"
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>

void investigate_scaling() {
    std::ofstream txt_file("scaling_results.txt");
    std::ofstream csv_file("scaling_results.csv");

    if (!txt_file.is_open() || !csv_file.is_open()) {
        std::cerr << "Error opening file for writing results.\n";
        return;
    }

    txt_file << "Investigating Scaling of the Method:\n";

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 10.0);

    for (int size = 1; size <= 30; ++size) {
        bool converged = false;
        int max_retries = 4;
        int retries = 0;

        while (!converged && retries < max_retries) {
            matrix H(size, size);
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j <= i; ++j) {
                    double value = dis(gen);
                    H.set(i, j, value);
                    H.set(j, i, value);
                }
            }

            auto start = std::chrono::high_resolution_clock::now();
            try {
                auto eigenpair = FindEigenPair(H, 1e-6, 2.0);
                auto end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> duration = end - start;

                txt_file << "Matrix size: " << size << "x" << size << ", Time: " << duration.count() << " seconds\n";
                csv_file << size << " " << duration.count() << "\n";

                converged = true;
            } catch (const std::exception &e) {
                txt_file << "Retrying due to non-convergence...\n";
                ++retries;
            }
        }

        if (!converged) {
            txt_file << "Stopping at size " << size << " due to repeated non-convergence.\n";
            throw std::runtime_error("Repeated non-convergence encountered.");
        }
    }

    txt_file.close();
    csv_file.close();
}
