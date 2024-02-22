#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>

// Define error function
double erf(double x) {
    if (x < 0) return -erf(-x);
    double a[] = {0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429};
    double t = 1.0 / (1.0 + 0.3275911 * x);
    double sum = t * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4]))));
    return 1.0 - sum * exp(-x * x);
}

// Define gamma function
double gamma(double x) {
    if (x < 0) return M_PI / sin(M_PI * x) / gamma(1 - x);
    if (x < 9) return gamma(x + 1) / x;
    double lngamma = x * log(x + 1 / (12 * x - 1 / (10 * x))) - x + log(2 * M_PI / x) / 2;
    return exp(lngamma);
}

// Define log gamma function
double lngamma(double x) {
    if (x <= 0) throw std::domain_error("lngamma: x <= 0");
    if (x < 9) return lngamma(x + 1) - log(x);
    return x * log(x + 1 / (12 * x - 1 / (10 * x))) - x + log(2 * M_PI / x) / 2;
}

// Main function
int main() {

    // Open files for writing
    std::ofstream gammaFile("gamma_tab.txt"), lngammaFile("lngamma_tab.txt"), erfFile("erf_data.txt"),
                  gammaDataFile("gamma_data.txt"), lngammaDataFile("lngamma_data.txt");

    if (!gammaFile.is_open() || !lngammaFile.is_open() || !erfFile.is_open() ||
        !gammaDataFile.is_open() || !lngammaDataFile.is_open()) {
        std::cerr << "Error opening file(s) for writing." << std::endl;
        return 1;
    }

    // Generate and write continuous values for erf
    for (double x = -3.0; x <= 3.0; x += 0.1) {
        double erfVal = erf(x);
        erfFile << x << " " << erfVal << std::endl;
    }

    // Generate and write continuous values for gamma
    for (double x = -10.0; x <= 10.0; x += 0.01) {
        double gammaVal = gamma(x);
        gammaDataFile << x << " " << gammaVal << std::endl;
    }

    // Generate and write continuous values for lngamma
    for (double x = 0.1; x <= 10.0; x += 0.01) {
        double lngammaVal;

        try {
            lngammaVal = lngamma(x);
        } catch (const std::exception& e) {
            std::cerr << "Exception for lngamma at x = " << x << ": " << e.what() << std::endl;
            continue;
        }

        lngammaDataFile << x << " " << lngammaVal << std::endl;
    }

    // Generate and write tabulated values for gamma and lngamma
    for (int n = 1; n <= 10; ++n) {
        double gammaVal = std::tgamma(n);
        double lngammaVal = std::lgamma(n);

        gammaFile << n << " " << gammaVal << std::endl;
        lngammaFile << n << " " << lngammaVal << std::endl;
    }

    gammaFile.close();
    lngammaFile.close();
    erfFile.close();
    gammaDataFile.close();
    lngammaDataFile.close();

    return 0;
}
