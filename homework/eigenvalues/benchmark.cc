#include "JDWCS.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>
#include <random>

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

void benchmarkJDWCS(int startSize, int endSize, std::ofstream& outFile) {
    for (int m = startSize; m <= endSize; ++m) {
        matrix A = rnd_symmetric_matrix(m);
        
        auto start = std::chrono::high_resolution_clock::now();
        JDWCS::cyclic(A); // Assuming this function performs the diagonalization
        auto end = std::chrono::high_resolution_clock::now();
        
        std::chrono::duration<double> elapsed = end - start;
        outFile << m << " " << elapsed.count() << "\n";
    }
}

int main() {
    const int nthreads = 4; // Number of threads to use
    const int nmax = 100; // Max matrix size
    const int nmin = 3; // Min matrix size
    const int step = (nmax - nmin + 1) / nthreads; // Calculate range of sizes for each thread

    std::vector<std::thread> threads;
    std::vector<std::ofstream> outFiles(nthreads);

    // Create threads and assign each a subrange of matrix sizes
    for (int i = 0; i < nthreads; ++i) {
        int start = nmin + i * step;
        int end = (i == nthreads - 1) ? nmax : start + step - 1; // Last thread gets the remainder
        outFiles[i].open("benchmark_" + std::to_string(i) + ".csv");
        threads.emplace_back(benchmarkJDWCS, start, end, std::ref(outFiles[i]));
    }

    // Join threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Close files
    for (auto& file : outFiles) {
        file.close();
    }

    // Aggregate results into a single file (Optional)
    std::ofstream finalFile("benchmark.csv");
    finalFile << "m, time\n";
    for (int i = 0; i < nthreads; ++i) {
        std::ifstream inFile("benchmark_" + std::to_string(i) + ".csv");
        std::string line;
        while (std::getline(inFile, line)) {
            finalFile << line << "\n";
        }
        inFile.close();
    }
    finalFile.close();

    return 0;
}