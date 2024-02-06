#include <iostream>
#include "vec.h"

int main() {
    // Create some vec objects
    vec v1(1.0, 2.0, 3.0);
    vec v2(4.0, 5.0, 6.0);

    // Perform operations and display results
    vec v3 = v1 + v2;
    vec v4 = v1 * 3.0;
    vec v5 = 2.0 * v2;
    vec v6 = v1 - v2;
    vec v7 = -v1;

    std::cout << "v1: " << v1.to_string() << std::endl;
    std::cout << "v2: " << v2.to_string() << std::endl;
    std::cout << "v3 = v1 + v2: " << v3.to_string() << std::endl;
    std::cout << "v4 = v1 * 3.0: " << v4.to_string() << std::endl;
    std::cout << "v5 = 2.0 * v2: " << v5.to_string() << std::endl;
    std::cout << "v6 = v1 - v2: " << v6.to_string() << std::endl;
    std::cout << "v7 = -v1: " << v7.to_string() << std::endl;

    double dotProduct = vec::dot(v1, v2);
    vec crossProduct = v1.cross(v2);
    double normV1 = v1.norm();

    std::cout << "Dot product (v1 dot v2): " << dotProduct << std::endl;
    std::cout << "Cross product (v1 cross v2): " << crossProduct.to_string() << std::endl;
    std::cout << "Norm of v1: " << normV1 << std::endl;

    // Test approximation
    bool approxCheck = vec::approx(v1, v2);
    std::cout << "Approximation check (v1 approx v2): " << (approxCheck ? "true" : "false") << std::endl;

    return 0;
}
