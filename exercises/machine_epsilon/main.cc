#include <iostream>
#include <limits> // Include for std::numeric_limits
#include <iomanip> // Include for std::scientific (for scientific notation)
#include "math_utils.h"

int main() {

    // Call the functions to find the max and min int
    int max_int = find_max_int();
    int min_int = find_min_int();

    // Call the functions to find the machine epsilon for double and float
    double epsilon_d = double_epsilon();
    float epsilon_f = float_epsilon();

    // Call the functions to find the sum of the series using method A and B
    int n = 1e6;
    double sumA = sum_method_A(n, epsilon_d);
    double sumB = sum_method_B(n, epsilon_d);

    // Define the same number in two different ways
    double d1 = 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1;
    double d2 = 8 * 0.1;

    // Compare the two numbers with a given accuracy
    bool approx_equal = approx(d1, d2);

    // Print the results and compare them to the expected values
    std::cout << "My max int = " << max_int << ", Should be = " << std::numeric_limits<int>::max() << std::endl;
    std::cout << "My min int = " << min_int << ", Should be = " << std::numeric_limits<int>::min() << std::endl;
    std::cout << "My double epsilon = " << epsilon_d << ", Should be = " << std::numeric_limits<double>::epsilon() << std::endl;
    std::cout << "My float epsilon = " << epsilon_f << ", Should be = " << std::numeric_limits<float>::epsilon() << std::endl;
    std::cout << "My sumA - 1 = " << std::scientific << sumA - 1 << ", Should be = " << std::scientific << n * (epsilon_d / 2.0) << std::endl;
    std::cout << "My sumB - 1 = " << std::scientific << sumB - 1 << ", Should be = " << std::scientific << n * (epsilon_d / 2.0) << std::endl;
    std::cout << "d1 == d2 ? => " << (d1 == d2) << ", Should be = 1" << std::endl;
    std::cout << "\nThe difference between sumA-1 and sumB-1 arises from floating-point arithmetic's limitations.\n"
             "In sumA, adding 1 before tiny values leads to precision loss due to rounding errors as tiny's effect\n"
             "becomes negligible. In sumB, adding tiny values before 1 accumulates more accurately. This demonstrates\n"
             "floating-point addition's non-associativity, where the order of operations affects the outcome due to\n"
             "limited precision.\n" << std::endl;
    std::cout << "approx(d1, d2) ? => " << approx_equal << ", Should be = 1" << std::endl;

    return 0;
}