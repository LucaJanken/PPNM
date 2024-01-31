#include <iostream>
#include <cmath>

int main() {
  double result = sqrt(2);
  double result2 = pow(2, 1.0/5.0);
  double result3 = exp(M_PI);
  double result4 = pow(M_PI, exp(1));
  std::cout << "sqrt(2) = " << result << " (should be approximately 1.4142)" << std::endl;
  std::cout << "2^(1/5) = " << result2 << " (should be approximately 1.1487)" << std::endl;
  std::cout << "exp(pi) = " << result3 << " (should be approximately 23.1407)" << std::endl;
  std::cout << "pi^(exp(1)) = " << result4 << " (should be approximately 22.4592)" << std::endl;
  return 0;
}
