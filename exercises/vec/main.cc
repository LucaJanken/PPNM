#include "vec.h"
#include <iostream>
#include <cmath>
#include <cassert>

void test_addition() {
    vec v1(1, 2, 3);
    vec v2(4, 5, 6);
    vec sum = v1 + v2;
    std::cout << "Test Addition: " << sum.to_string() << " == 5 7 9" << std::endl;
    assert(sum.approx(vec(5, 7, 9)));
}

void test_subtraction() {
    vec v1(10, 20, 30);
    vec v2(5, 5, 5);
    vec diff = v1 - v2;
    std::cout << "Test Subtraction: " << diff.to_string() << " == 5 15 25" << std::endl;
    assert(diff.approx(vec(5, 15, 25)));
}

void test_scalar_multiplication_left() {
    double scalar = 2.0;
    vec v(1, 2, 3);
    vec result = scalar * v; // Scalar multiplication from the left
    std::cout << "Test Scalar Multiplication (Left): " << result.to_string() << " == 2 4 6" << std::endl;
    assert(result.approx(vec(2, 4, 6)));
}

void test_scalar_multiplication_right() {
    vec v(1, 2, 3);
    double scalar = 2.0;
    vec result = v * scalar; // Scalar multiplication from the right
    std::cout << "Test Scalar Multiplication (Right): " << result.to_string() << " == 2 4 6" << std::endl;
    assert(result.approx(vec(2, 4, 6)));
}

void test_dot_product() {
    vec v1(1, 2, 3);
    vec v2(4, 5, 6);
    double dot = v1.dot(v2);
    std::cout << "Test Dot Product: " << dot << " == 32" << std::endl;
    assert(vec::approx(dot, 32));
}

void test_cross_product() {
    vec v1(1, 2, 3);
    vec v2(4, 5, 6);
    vec cross = v1.cross(v2);
    std::cout << "Test Cross Product: " << cross.to_string() << " == -3 6 -3" << std::endl;
    assert(cross.approx(vec(-3, 6, -3)));
}

void test_norm() {
    vec v(3, 4, 0);
    double norm = v.norm();
    std::cout << "Test Norm: " << norm << " == 5" << std::endl;
    assert(vec::approx(norm, 5));
}

void test_static_dot_product() {
    vec v1(1, 2, 3);
    vec v2(4, 5, 6);
    double dot_static = vec::dot(v1, v2);
    std::cout << "Test Static Dot Product: " << dot_static << " == 32" << std::endl;
    assert(vec::approx(dot_static, 32));
}

void test_static_cross_product() {
    vec v1(1, 2, 3);
    vec v2(4, 5, 6);
    vec cross_static = vec::cross(v1, v2);
    std::cout << "Test Static Cross Product: " << cross_static.to_string() << " == -3 6 -3" << std::endl;
    assert(cross_static.approx(vec(-3, 6, -3)));
}

void test_static_norm() {
    vec v(3, 4, 0);
    double norm_static = vec::norm(v);
    std::cout << "Test Static Norm: " << norm_static << " == 5" << std::endl;
    assert(vec::approx(norm_static, 5));
}

void test_static_approximation() {
    vec v1(1 + 1e-10, 2 + 1e-10, 3 + 1e-10);
    vec v2(1, 2, 3);
    bool is_approx_equal = vec::approx(v1, v2);
    std::cout << "Test Static Approximation: " << std::boolalpha << is_approx_equal << " == true" << std::endl;
    assert(is_approx_equal);
}

void run_tests() {
    // Instance methods
    test_addition();
    test_subtraction();
    test_scalar_multiplication_left();
    test_scalar_multiplication_right();
    test_dot_product();
    test_cross_product();
    test_norm();
    // Static methods
    test_static_dot_product();
    test_static_cross_product();
    test_static_norm();
    test_static_approximation();
    std::cout << "All tests passed successfully." << std::endl;
}

int main() {
    run_tests();
    return 0;
}
