#include "vec.h"
#include <cmath> // For sqrt and abs functions

// Default constructor
vec::vec() : x(0), y(0), z(0) {}

// Parameterized constructor
vec::vec(double x, double y, double z) : x(x), y(y), z(z) {}

// Overload the addition operator to add two vectors
vec vec::operator+(const vec& v) const {
    return vec(x + v.x, y + v.y, z + v.z);
}

// Overload the subtraction operator to subtract two vectors
vec vec::operator-(const vec& v) const {
    return vec(x - v.x, y - v.y, z - v.z);
}

// Overload the unary negation operator to negate a vector
vec vec::operator-() const {
    return vec(-x, -y, -z);
}

// Overload the multiplication operator to scale a vector by a scalar
vec vec::operator*(double c) const {
    return vec(x * c, y * c, z * c);
}

// Friend function to allow scalar multiplication on the left
vec operator*(double c, const vec& v) {
    return v * c;
}

// Dot product
double vec::dot(const vec& other) const {
    return x * other.x + y * other.y + z * other.z;
}

// Static method for dot product
double vec::dot(const vec& v, const vec& w) {
    return v.dot(w);
}

// Cross product
vec vec::cross(const vec& other) const {
    return vec(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
}

// Static method for cross product
vec vec::cross(const vec& v, const vec& w) {
    return v.cross(w);
}

// Norm of the vector
double vec::norm() const {
    return std::sqrt(dot(*this));
}

// Static method for norm
double vec::norm(const vec& v) {
    return v.norm();
}

// Static approximation method for two doubles
bool vec::approx(double a, double b, double acc, double eps) {
    if (std::abs(a - b) < acc) {
        return true;
    }
    if (std::abs(a - b) < (std::abs(a) + std::abs(b)) * eps) {
        return true;
    }
    return false;
}

// Approximation method for two vec objects
bool vec::approx(const vec& other, double acc, double eps) const {
    if (!approx(x, other.x, acc, eps)) {
        return false;
    }
    if (!approx(y, other.y, acc, eps)) {
        return false;
    }
    if (!approx(z, other.z, acc, eps)) {
        return false;
    }
    return true;
}

// Static approximation method for two vec objects
bool vec::approx(const vec& u, const vec& v, double acc, double eps) {
    return u.approx(v, acc, eps);
}

// Print method for debugging
void vec::print(const std::string& s) const {
    std::cout << s << x << " " << y << " " << z << std::endl;
}

// Override to_string method
std::string vec::to_string() const {
    return std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z);
}
