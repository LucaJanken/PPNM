#include <iostream>
#include <thread>
#include <vector>
#include <string>
#include <numeric> 
#include <cstdlib> 

#ifdef USE_OPENMP
#include <omp.h>
#endif

struct Data {
    int a, b;
    double sum;
    Data() : a(0), b(0), sum(0) {}
};

void harm(Data* data) {
    data->sum = 0;
    for (int i = data->a; i < data->b; ++i) {
        data->sum += 1.0 / i;
    }
}

int main(int argc, char* argv[]) {
    int nthreads = 1;
    int nterms = 1e8; // Default values

    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto pos = arg.find(':');
        if (pos != std::string::npos) {
            std::string key = arg.substr(0, pos);
            std::string value = arg.substr(pos + 1);
            if (key == "-threads") {
                nthreads = std::atoi(value.c_str());
            } else if (key == "-terms") {
                nterms = static_cast<int>(std::stod(value));
            }
        }
    }

#ifdef USE_OPENMP
double totalSum = 0.0;
    #ifdef USE_TLS
    // Correct approach using reduction - simulates TLS by ensuring thread-safe updates to totalSum
    #pragma omp parallel for reduction(+:totalSum) num_threads(nthreads)
    for (int i = 1; i <= nterms; ++i) {
        totalSum += 1.0 / i;
    }
    #else
    // Intentionally incorrect approach without reduction - likely to cause race conditions
    #pragma omp parallel for num_threads(nthreads)
    for (int i = 1; i <= nterms; ++i) {
        totalSum += 1.0 / i;
    }
    #endif
#else
    // Prepare data objects
    std::vector<Data> params(nthreads);
    for (int i = 0; i < nthreads; ++i) {
        params[i].a = 1 + nterms / nthreads * i;
        params[i].b = 1 + nterms / nthreads * (i + 1);
    }
    params.back().b = nterms + 1; // Adjust the endpoint

    // Create and run threads
    std::vector<std::thread> threads;
    for (int i = 0; i < nthreads; ++i) {
        threads.emplace_back(harm, &params[i]);
    }

    // Join the threads
    for (auto& thread : threads) {
        thread.join();
    }

    // Calculate the total sum
    double totalSum = 0;
    for (const auto& p : params) {
        totalSum += p.sum;
    }
#endif

    std::cout << "Total Harmonic Sum: " << totalSum << std::endl;

    return 0;
}
