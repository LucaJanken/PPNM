#include "quadrature.h"

// Adjusted integrate signature
std::tuple<double, double, int> integrate(const std::function<double(double)>& f, double a, double b, double delta, double eps, double f2, double f3) {
    int evaluations = 0; // Initialize the evaluations count
    double h = b - a;
    if (std::isnan(f2)) {
        f2 = f(a + 2 * h / 6); evaluations++;
        f3 = f(a + 4 * h / 6); evaluations++;
    }
    double f1 = f(a + h / 6); evaluations++;
    double f4 = f(a + 5 * h / 6); evaluations++;
    double Q = (2 * f1 + f2 + f3 + 2 * f4) / 6 * h;
    double q = (f1 + f2 + f3 + f4) / 4 * h;
    double err = std::abs(Q - q);

    if (err <= delta + eps * std::abs(Q)) {
        return std::make_tuple(Q, err, evaluations);
    } else {
        auto left = integrate(f, a, (a + b) / 2, delta / std::sqrt(2), eps, f1, f2);
        auto right = integrate(f, (a + b) / 2, b, delta / std::sqrt(2), eps, f3, f4);
        double final_result = std::get<0>(left) + std::get<0>(right);
        double final_error = std::get<1>(left) + std::get<1>(right);
        int total_evaluations = evaluations + std::get<2>(left) + std::get<2>(right); // Accumulate evaluations from both sides and current
        return std::make_tuple(final_result, final_error, total_evaluations);
    }
}

std::tuple<double, double, int> transint(const std::function<double(double)>& f, double a, double b, double delta, double eps) {
    // Initial evaluations count for transint
    int evaluations = 0;
    std::function<double(double)> ft = [&f, a, b, &evaluations](double theta) -> double {
        evaluations++; // Increment evaluations here because each call to ft counts
        return f((a + b) / 2.0 + cos(theta) * (b - a) / 2.0) * sin(theta) * (b - a) / 2.0;
    };

    // Integrate the transformed function over [0, π], without initial evaluations since it's now internal to ft
    return integrate(ft, 0, M_PI, delta, eps, evaluations);
}

// Transformation for handling infinite limits
double transform(double x) {
    return x / (1.0 - std::abs(x));
}
// Adjusted integrate_infinite signature
std::tuple<double, double, int> integrate_infinite(const std::function<double(double)>& f, double a, double b, double delta, double eps) {
    // Check for infinite limits and apply transformations if necessary
    std::function<double(double)> transformed_f;
    double transformed_a = a, transformed_b = b;

    if (std::isinf(a) || std::isinf(b)) {
        transformed_f = [f](double x) {
            double t = transform(x); // Apply transformation
            double dt_dx = 1.0 / ((1.0 - std::abs(x)) * (1.0 - std::abs(x))); // Derivative of the transformation
            return f(t) * dt_dx; // Adjusted integrand
        };
        transformed_a = (a == -std::numeric_limits<double>::infinity()) ? -1.0 : a; // Transform to -1 if a is negative infinity
        transformed_b = (b == std::numeric_limits<double>::infinity()) ? 1.0 : b; // Transform to 1 if b is positive infinity
    } else {
        transformed_f = f; // Use the original function if limits are finite
    }

    // Use the existing integrate function for the actual computation
    return integrate(transformed_f, transformed_a, transformed_b, delta, eps);
}