#ifndef HAVE_RKODE_H
#define HAVE_RKODE_H
#include "../matrix/vector.h"
#include <vector>
#include <functional>
#include <cmath>

// Correct the function name to match the implementation if necessary
std::pair<vector, vector> rkstep4 (
    const std::function<vector(double, const vector&)>& f,  // Function dy/xy = f(x, y)
    double x,                                               // Current x
    const vector& y,                                        // Current y(x)
    double h                                                // Step size
);

// Corrected driver method signature
std::pair<vector, std::vector<vector>> driver(
    const std::function<vector(double, const vector&)>& F, // the function f from dy/dx=f(x,y)
    std::pair<double, double> interval, // (start-point, end-point)
    vector ystart,               // y(start-point)
    double h = 0.01,            // initial step-size
    double acc = 0.01,           // absolute accuracy goal
    double eps = 0.01,           // relative accuracy goal
    int nmax = 9999               // maximum number of steps
);

// Interpolation driver method
std::function<vector(double)> driver_interp(
    const std::function<vector(double, const vector&)>& F,               // the function f from dy/dx=f(x,y)
    std::pair<double, double> interval, // (start-point, end-point)
    vector ystart,               // y(start-point)
    double h = 0.01,            // initial step-size
    double acc = 0.01,           // absolute accuracy goal
    double eps = 0.01,           // relative accuracy goal
    int nmax = 9999               // maximum number of steps (added to match the implementation)
);


#endif
