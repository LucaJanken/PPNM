#include <stdexcept>
#include <cmath>
#include "../matrix/vector.h"
#include "spline.h"

// Implement the static binsearch method from Spline
int Spline::binsearch(const vector& x, double z) {
   /* locates the interval for z by bisection */
    if (z < x[0] || z > x[x.size-1]) {
        throw std::invalid_argument("z must be within the range of x");
    }
    int i = 0;
    int j = x.size - 1;
    while (j - i > 1) {
        int mid = (i + j) / 2;
        if (z > x[mid]) {
            i = mid;
        } else {
            j = mid;
        }
    }
    return i;
}

// LSpline Constructor
LSpline::LSpline(const vector& xs, const vector& ys) : x(xs), y(ys) {}

// LSpline evaluate method
double LSpline::evaluate(double z) const {
    int i = binsearch(x, z);
    double dx = x[i+1] - x[i];
    if (!(dx > 0)) { // Throw exception
        throw std::invalid_argument("x must be strictly increasing");
    }
    double dy = y[i+1] - y[i];
    return y[i] + dy * (z - x[i]) / dx;
}

// LSpline integral method
double LSpline::integral(double z) const {
    int i = binsearch(x, z);
    double sum = 0;
    for (int j = 0; j < i; j++) {
        sum += (x[j+1]-x[j]) * (y[j+1]+y[j]) / 2;
    }
    sum += (z-x[i]) * (y[i]+y[i+1]) / 2;
    return sum;
}


// QSpline Constructor
QSpline::QSpline(const vector& xs, const vector& ys) {
    if (xs.size != ys.size) {
        throw std::runtime_error("X and Y vectors must have the same size");
    }

    int n = xs.size;
    this->x = xs;
    this->y = ys;

    vector b_temp(n - 1), c_temp(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        b_temp[i] = 0;
        c_temp[i] = 0;
    }
    this->b = b_temp;
    this->c = c_temp;
    
    vector p(n - 1), h(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        h[i] = this->x[i + 1] - this->x[i];
        p[i] = (this->y[i + 1] - this->y[i]) / h[i];
    }

    // Forward recursion for c's
    this->c[0] = 0.0;
    for (int i = 0; i < n - 2; ++i) {
        this->c[i + 1] = (p[i + 1] - p[i] - this->c[i] * h[i]) / h[i + 1];
    }
    this->c[n - 2] /= 2.0;

    // Backward recursion for c's, adjust c based on forward recursion results
    for (int i = n - 3; i >= 0; --i) {
        this->c[i] = (p[i + 1] - p[i] - this->c[i + 1] * h[i + 1]) / h[i];
    }

    // Calculate b's using c's
    for (int i = 0; i < n - 1; ++i) {
        this->b[i] = p[i] - this->c[i] * h[i];
    }
}

// QSpline evaluate method
double QSpline::evaluate(double z) const {
    int i = binsearch(x, z);
    double dx = z - x[i];
    return y[i] + b[i] * dx + c[i] * dx * dx;
}

// QSpline derivative method
double QSpline::derivative(double z) const {
    (void) z; // Suppress unused parameter warning
    int i = binsearch(x, z);
    double dx = z - x[i];
    return b[i] + 2 * c[i] * dx;
}

// QSpline integral method
double QSpline::integral(double z) const {
    int i = binsearch(x, z);
    double sum = 0;
    for (int j = 0; j < i; j++) {
        double dx = x[j+1] - x[j];
        sum += y[j] * dx + b[j] * dx * dx / 2 + c[j] * dx * dx * dx / 3;
    }
    double dx = z - x[i];
    sum += y[i] * dx + b[i] * dx * dx / 2 + c[i] * dx * dx * dx / 3;
    return sum;
}

// CSpline Constructor
CSpline::CSpline(const vector& xs, const vector& ys) {
    if (xs.size != ys.size) {
        throw std::runtime_error("X and Y vectors must have the same size");
    }
    int n;
    n = xs.size;
    x = xs;
    y = ys;

    b.resize(n);
    c.resize(n, 0.0);
    d.resize(n-1, 0.0);
    vector h(n-1), p(n-1);

    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
        if (h[i] <= 0) {
            throw std::runtime_error("x values must be strictly increasing");
        }
        p[i] = (y[i + 1] - y[i]) / h[i];
    }

    // Initializing vectors for tridiagonal solver
    vector alpha(n), l(n), mu(n), z(n);
    alpha[0] = 0.0;
    l[0] = 2.0;
    mu[0] = 0.0;
    z[0] = 0.0;

    for (int i = 1; i < n - 1; ++i) {
        alpha[i] = 3.0 * (p[i] - p[i - 1]);
        l[i] = 4.0 - mu[i - 1];
        mu[i] = 1.0 / l[i];
        z[i] = (alpha[i] - z[i - 1]) / l[i];
    }

    l[n - 1] = 1.0;
    z[n - 1] = 0.0;
    c[n - 1] = 0.0;

    // Backward pass for c
    for (int i = n - 2; i >= 0; --i) {
        c[i] = z[i] - mu[i] * c[i + 1];
    }

    // Solving for b and d
    for (int i = 0; i < n - 1; ++i) {
        b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0;
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
    }
}

// CSpline evaluate method
double CSpline::evaluate(double z) const {
    int i = binsearch(x, z);
    double dx = z - x[i];
    return y[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
}

// CSpline derivative method
double CSpline::derivative(double z) const {
    (void) z; // Suppress unused parameter warning
    int i = binsearch(x, z);
    double dx = z - x[i];
    return b[i] + 2 * c[i] * dx + 3 * d[i] * dx * dx;
}

// CSpline integral method
double CSpline::integral(double z) const {
    int i = binsearch(x, z);
    double sum = 0;
    for (int j = 0; j < i; j++) {
        double dx = x[j+1] - x[j];
        sum += y[j] * dx + b[j] * dx * dx / 2 + c[j] * dx * dx * dx / 3 + d[j] * dx * dx * dx * dx / 4;
    }
    double dx = z - x[i];
    sum += y[i] * dx + b[i] * dx * dx / 2 + c[i] * dx * dx * dx / 3 + d[i] * dx * dx * dx * dx / 4;
    return sum;
}
