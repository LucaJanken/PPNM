#include <iostream>
#include <cmath>
#include <complex>
#include "vec.h"

int main() {
    // Expressions
    std::complex<double> i(0, 1);
    auto sqrt_neg1 = std::sqrt(std::complex<double>(-1, 0));
    auto sqrt_i = std::sqrt(i);
    auto e_i = std::exp(i);
    auto e_i_pi = std::exp(i * M_PI);
    auto i_i = std::pow(i, i);
    auto ln_i = std::log(i);
    auto sin_i_pi = std::sin(i * M_PI);
    auto sinh_i_pi = std::sinh(i * M_PI); 
    auto cosh_i_pi = std::cosh(i * M_PI); 

    // Manually calculated results
    std::complex<double> manual_sqrt_neg1(0, 1); // √-1 = i
    std::complex<double> manual_sqrt_i(std::sqrt(0.5), std::sqrt(0.5)); // √i = (1/√2) + i(1/√2)
    std::complex<double> manual_e_i(std::cos(1), std::sin(1)); // e^i = cos(1) + i*sin(1)
    std::complex<double> manual_e_i_pi(-1, 0); // e^(i*pi) = -1
    std::complex<double> manual_i_i(std::exp(-M_PI / 2), 0); // i^i = e^(-π/2)
    std::complex<double> manual_ln_i(0, M_PI / 2); // ln(i) = iπ/2
    std::complex<double> manual_sin_i_pi(0, std::sinh(M_PI)); // sin(i*pi) = i*sinh(pi)
    std::complex<double> expected_sinh_i_pi(0, 0); // sinh(i*pi) = 0
    std::complex<double> expected_cosh_i_pi(-1, 0); // cosh(i*pi) = -1

    // Output results
    std::cout << "sqrt(-1) = " << sqrt_neg1 << " (expected: " << manual_sqrt_neg1 << ")" << std::endl;
    std::cout << "sqrt(i) = " << sqrt_i << " (expected: " << manual_sqrt_i << ")" << std::endl;
    std::cout << "exp(i) = " << e_i << " (expected: " << manual_e_i << ")" << std::endl;
    std::cout << "exp(i * pi) = " << e_i_pi << " (expected: " << manual_e_i_pi << ")" << std::endl;
    std::cout << "i^i = " << i_i << " (expected: " << manual_i_i << ")" << std::endl;
    std::cout << "ln(i) = " << ln_i << " (expected: " << manual_ln_i << ")" << std::endl;
    std::cout << "sin(i * pi) = " << sin_i_pi << " (expected: " << manual_sin_i_pi << ")" << std::endl;
    std::cout << "sinh(i * pi) = " << sinh_i_pi << " (expected: " << expected_sinh_i_pi << ")" << std::endl;
    std::cout << "cosh(i * pi) = " << cosh_i_pi << " (expected: " << expected_cosh_i_pi << ")" << std::endl;
    

    // Comparison function for complex numbers using approx
    auto compare_complex = [](const std::complex<double>& computed, const std::complex<double>& expected) {
        return vec::approx(computed.real(), expected.real()) && vec::approx(computed.imag(), expected.imag());
    };

    // Perform comparisons
    std::cout << "Comparing computed results with expected values:" << std::endl;
    std::cout << "sqrt(-1) comparison: " << (compare_complex(sqrt_neg1, manual_sqrt_neg1) ? "PASS" : "FAIL") << std::endl;
    std::cout << "sqrt(i) comparison: " << (compare_complex(sqrt_i, manual_sqrt_i) ? "PASS" : "FAIL") << std::endl;
    std::cout << "exp(i) comparison: " << (compare_complex(e_i, manual_e_i) ? "PASS" : "FAIL") << std::endl;
    std::cout << "exp(i * pi) comparison: " << (compare_complex(e_i_pi, manual_e_i_pi) ? "PASS" : "FAIL") << std::endl;
    std::cout << "i^i comparison: " << (compare_complex(i_i, manual_i_i) ? "PASS" : "FAIL") << std::endl;
    std::cout << "ln(i) comparison: " << (compare_complex(ln_i, manual_ln_i) ? "PASS" : "FAIL") << std::endl;
    std::cout << "sin(i * pi) comparison: " << (compare_complex(sin_i_pi, manual_sin_i_pi) ? "PASS" : "FAIL") << std::endl;
    std::cout << "sinh(i * pi) comparison: " << (compare_complex(sinh_i_pi, expected_sinh_i_pi) ? "PASS" : "FAIL") << std::endl;
    std::cout << "cosh(i * pi) comparison: " << (compare_complex(cosh_i_pi, expected_cosh_i_pi) ? "PASS" : "FAIL") << std::endl;

    return 0;
}