#include <iostream>
#include <cmath>
#include "sfuns.h"

int main() {
  double result = sqrt(2);
  double result2 = pow(2, 1.0/5.0);
  double result3 = exp(M_PI);
  double result4 = pow(M_PI, exp(1));

  std::cout << "sqrt(2) = " << result << " (should be approximately 1.4142)" << std::endl;
  std::cout << "2^(1/5) = " << result2 << " (should be approximately 1.1487)" << std::endl;
  std::cout << "exp(pi) = " << result3 << " (should be approximately 23.1407)" << std::endl;
  std::cout << "pi^(exp(1)) = " << result4 << " (should be approximately 22.4592)" << std::endl;

  std::cout << "Gamma(1) = " << my_gamma(1) << " (should be exactly 1)" << std::endl;
  std::cout << "Gamma(2) = " << my_gamma(2) << " (should be exactly 1)" << std::endl;
  std::cout << "Gamma(3) = " << my_gamma(3) << " (should be exactly 2)" << std::endl;
  std::cout << "Gamma(10) = " << my_gamma(10) << " (should be exactly 362880)" << std::endl;

  std::cout << "lnGamma(1) = " << my_lngamma(1) << " (should be exactly 0)" << std::endl;
  std::cout << "lnGamma(2) = " << my_lngamma(2) << " (should be exactly 0)" << std::endl;
  std::cout << "lnGamma(3) = " << my_lngamma(3) << " (should be approximately 0.6931)" << std::endl;
  std::cout << "lnGamma(10) = " << my_lngamma(10) << " (should be approximately 12.8018)" << std::endl;
  
  return 0;
}
