#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
#include <tuple>
#include <fstream>
#include "../matrix/vector.h"
#include "../matrix/matrix.h"
#include "lsfit.h"

// Util function to convert vector to log space
vector log_vector(const vector& v) {
    vector result(v.size);
    for(size_t i = 0; i < v.size; i++) {
        result.set(i, log(v.get(i)));
    }
    return result;
}

// Util function to calculate uncertainties
vector calc_uncertainties(const vector& v, const vector& dv) {
    vector result(v.size);
    for(size_t i = 0; i < v.size; i++) {
        result.set(i, dv.get(i) / v.get(i));
    }
    return result;
}

int main() {
    // Declare vectors with sizes
    vector t(9), y(9), dy(9);

    // Read data from file
    std::ifstream file("datafile.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open datafile.csv" << std::endl;
        return -1;
    }
    int index = 0;
    double time, activity, uncertainty;
    while (file >> time >> activity >> uncertainty) {
        t.set(index, time);
        y.set(index, activity);
        dy.set(index, uncertainty);
        index++;
    }
    file.close();

    // Convert to log space
    vector lny = log_vector(y);
    vector dlny = calc_uncertainties(y, dy);

    // Define basis functions
    std::vector<std::function<double(double)>> fs = {
        [](double) -> double { return 1; },
        [](double x) -> double { return -x; }
    };

    // Perform least squares fit
    vector c;
    matrix cov;
    std::tie(c, cov) = lsfit(fs, t, lny, dlny); 

    // Print statements
    c.print("c:");
    std::cout << "Uncertainties: " << sqrt(cov.get(0, 0)) << " " << sqrt(cov.get(1, 1)) << std::endl;
    
    double tau = log(2) / c.get(1);
    double dtau = log(2) / (c.get(1) * c.get(1)) * sqrt(cov.get(1, 1));
    
    std::cout << "a = " << exp(c.get(0)) << std::endl;
    std::cout << "tau = " << tau << " pm " << dtau << " days, table value: 3.6313(14) days" << std::endl;
    
    cov.print("covariance matrix:");

    return 0;
}
