#include <iostream>
#include <cmath>
#include "sfuns.h"

int main() {

  // Assigning the results of the mathematical functions to variables (using cmath library functions)
  double sqrt_of_two = sqrt(2);
  double fifth_root_of_two = pow(2, 1.0/5.0);
  double exp_of_pi = exp(M_PI);
  double pi_to_the_power_of_e = pow(M_PI, exp(1));

  // Printing the results and comparing them to the expected values
  std::cout << "sqrt(2) = " << sqrt_of_two << " (should be approximately 1.4142)" << std::endl;
  std::cout << "2^(1/5) = " << fifth_root_of_two << " (should be approximately 1.1487)" << std::endl;
  std::cout << "exp(pi) = " << exp_of_pi << " (should be approximately 23.1407)" << std::endl;
  std::cout << "pi^(exp(1)) = " << pi_to_the_power_of_e << " (should be approximately 22.4592)" << std::endl;
  
  // Printing the results of the Gamma function and comparing them to the expected values
  std::cout << "Gamma(1) = " << my_gamma(1) << " (should be exactly 1)" << std::endl;
  std::cout << "Gamma(2) = " << my_gamma(2) << " (should be exactly 1)" << std::endl;
  std::cout << "Gamma(3) = " << my_gamma(3) << " (should be exactly 2)" << std::endl;
  std::cout << "Gamma(10) = " << my_gamma(10) << " (should be exactly 362880)" << std::endl;

  // Printing the results of the lnGamma function and comparing them to the expected values
  std::cout << "lnGamma(1) = " << my_lngamma(1) << " (should be exactly 0)" << std::endl;
  std::cout << "lnGamma(2) = " << my_lngamma(2) << " (should be exactly 0)" << std::endl;
  std::cout << "lnGamma(3) = " << my_lngamma(3) << " (should be approximately 0.6931)" << std::endl;
  std::cout << "lnGamma(10) = " << my_lngamma(10) << " (should be approximately 12.8018)" << std::endl;
  
  return 0;
}
