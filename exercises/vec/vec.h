#ifndef VEC_H
#define VEC_H

#include <iostream>
#include <string>

class vec {
    public:
        static constexpr double DEFAULT_ACC = 1e-9;
        static constexpr double DEFAULT_EPS = 1e-9;

        double x, y, z; // Declare x, y, and z as components of the vector

        // Default constructor
        vec();

        // Parameterized constructor
        vec(double x, double y, double z);

        // Overload the addition operator to add two vectors
        vec operator+(const vec& v) const;

        // Overload the subtraction operator to subtract two vectors
        vec operator-(const vec& v) const;

        // Overload the unary negation operator to negate a vector
        vec operator-() const;

        // Overload the multiplication operator to scale a vector by a scalar
        vec operator*(double c) const;

        // Friend function to allow scalar multiplication on the left
        friend vec operator*(double c, const vec& v);

        // Dot product
        double dot(const vec& other) const;

        // Static method for dot product
        static double dot(const vec& v, const vec& w);

        // Cross product
        vec cross(const vec& other) const;

        // Static method for cross product
        static vec cross(const vec& v, const vec& w);

        // Norm of the vector
        double norm() const;

        // Static method for norm
        static double norm(const vec& v);

        // Approximation method (to compare two doubles for approximate equality)
        static bool approx(double a, double b, double acc = DEFAULT_ACC, double eps = DEFAULT_EPS);

        // Approximation method for two vec objects
        bool approx(const vec& other, double acc = DEFAULT_ACC, double eps = DEFAULT_EPS) const;

        // Static approximation method for two vec objects
        static bool approx(const vec& u, const vec& v, double acc = DEFAULT_ACC, double eps = DEFAULT_EPS);

        // Print method for debugging
        void print(const std::string& s = "") const;

        // Override to_string method
        std::string to_string() const;
};

#endif // VEC_H
