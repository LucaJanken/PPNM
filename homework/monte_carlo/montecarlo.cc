#include "montecarlo.h"

std::pair<double, double> plainmc(std::function<double(vector&)> f, const vector& a, const vector& b, int N) {
    int dim = a.size;
    double V = 1.0;
    for(int i = 0; i < dim; ++i) V *= (b[i] - a[i]);
    
    double sum = 0.0, sum2 = 0.0;
    vector x(dim);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for(int i = 0; i < N; ++i){
        for(int k = 0; k < dim; ++k) x[k] = a[k] + dis(gen)*(b[k] - a[k]);
        double fx = f(x);
        sum += fx;
        sum2 += fx*fx;
    }

    double mean = sum / N;
    double sigma = sqrt(sum2 / N - mean * mean);
    return std::make_pair(mean * V, sigma * V / sqrt(N));
}

double corput(int n , int b){
    double q = 0, bk = (double) 1 / b;
    while(n > 0) { 
        q += (n % b) * bk ; n /= b ; bk /= b; 
    }
    return q; 
}

void halton(int n, int d, double *x) {
    static int base[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};
    static int maxd = sizeof(base) / sizeof(int);
    assert(d <= maxd);
    for (int i = 0; i < d; i++) x[i] = corput(n, base[i]);
}

void lattice_rule(int n, int d, double *x) {
    double primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};
    int numPrimes = sizeof(primes) / sizeof(primes[0]);
    for (int i = 0; i < numPrimes; i++) {
        primes[i] = sqrt(primes[i]);
    }
    if (d > numPrimes) {
        std::cerr << "Error: d > numPrimes" << std::endl;
        exit(1);
    }
    for (int i = 0; i < d; i++) {
        // The fractional part of (sqrt(prime[i]) * n)
        x[i] = fmod(primes[i] * n, 1.0);
    }
}

std::pair<double, double> quasi_mc(std::function<double(vector&)> f, const vector& a, const vector& b, int N) {
    int dim = a.size;
    double V = 1.0;
    for (int i = 0; i < dim; ++i) {
        V *= (b[i] - a[i]);
    }

    double sum1 = 0.0, sum2 = 0.0;
    vector x1(dim), x2(dim);
    
    for (int i = 0; i < N; ++i) {
        // First quasi-random sequence (e.g., Halton)
        halton(i, dim, &x1[0]);
        for (int k = 0; k < dim; ++k) x1[k] = a[k] + x1[k] * (b[k] - a[k]);
        sum1 += f(x1);

        // Second quasi-random sequence (e.g., lattice rule)
        lattice_rule(i, dim, &x2[0]);
        for (int k = 0; k < dim; ++k) x2[k] = a[k] + x2[k] * (b[k] - a[k]);
        sum2 += f(x2);
    }

    double mean1 = sum1 / N;
    double mean2 = sum2 / N;
    
    double error = std::abs(mean1 - mean2) * V;
    double area = (mean1 * V + mean2 * V) / 2; // Average the estimates for the final area

    return std::make_pair(area, error);
}

std::pair<double, double> stratified_mc(std::function<double(vector&)> f, const vector& a, const vector& b, int N, int nmin) {
    if (N <= nmin) {
        return plainmc(f, a, b, N); // Use plain Monte Carlo for small N
    }
    int dim = a.size;
    double V = 1.0;
    for (int i = 0; i < dim; ++i) {
        V *= (b[i] - a[i]);
    }
    double sum = 0.0, sum2 = 0.0;
    vector x(dim);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    matrix sums(2, dim), sums2(2, dim), vars(2, dim);
    for (int i = 0; i < nmin; ++i) {
        for (int k = 0; k < dim; ++k) {
            x[k] = a[k] + dis(gen) * (b[k] - a[k]);
        }
        double fx = f(x);
        sum += fx;
        sum2 += fx * fx;
        for (int k = 0; k < dim; ++k) {
            if(x[k] < (a[k] + b[k]) / 2) {
                sums(0, k) += fx;
                sums2(0, k) += fx * fx;
            } else {
                sums(1, k) += fx;
                sums2(1, k) += fx * fx;
            }
        }
    }
    double mean = sum / nmin;
    double sigma = V * sqrt(sum2 / nmin - mean * mean) / sqrt(nmin);
    int wdim = 0;
    double max_var = 0.0;
    for (int k = 0; k < dim; ++k) {
        vars(0, k) = (sums2(0, k) - sums(0, k) * sums(0, k) / nmin) / nmin;
        vars(1, k) = (sums2(1, k) - sums(1, k) * sums(1, k) / nmin) / nmin;
        if (vars(0, k) > max_var) {
            max_var = vars(0, k);
            wdim = k;
        }
        if (vars(1, k) > max_var) {
            max_var = vars(1, k);
            wdim = k;
        }
    }
    vector a1 = a, b1 = b;
    a1[wdim] = (a[wdim] + b[wdim]) / 2;
    b1[wdim] = (a[wdim] + b[wdim]) / 2;
    int N_up = 0;
    while (N_up + 1 < (N - nmin) / (1 + vars(0, wdim) / vars(1, wdim))) {
        N_up++;
    }
    int N_down = N - nmin - N_up;
    std::pair<double, double> result_down = stratified_mc(f, a, b1, N_down, nmin);
    std::pair<double, double> result_up = stratified_mc(f, a1, b, N_up, nmin);
    double total_integral = ((result_down.first + result_up.first) * (N - nmin) + V * mean * nmin) / N;
    double total_error = sqrt(pow(result_down.second * (N - nmin) / N, 2) + pow(result_up.second * (N - nmin) / N, 2) + pow(sigma * nmin / N, 2));
    return std::make_pair(total_integral, total_error); 
}