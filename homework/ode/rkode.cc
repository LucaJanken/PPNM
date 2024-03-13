#include "rkode.h"
#include <iostream>
#include "../splines/spline.h"

// RK12 (Euler/Midpoint) method
std::pair<vector, vector> rkstep4 (
    const std::function<vector(double, const vector&)>& f,  // Function dy/xy = f(x, y)
    double x,                                               // Current x
    const vector& y,                                        // Current y(x)
    double h                                                // Step size
) {
    vector k0 = f(x, y);                                    // Lower order step (Euler)
    vector k1 = f(x + h / 2, y + k0 * (h / 2));             // Higher order step (Midpoint)
    vector k2 = f(x + h / 2, y + k1 * (h / 2));
    vector k3 = f(x + h, y + k2 * h);
    vector k = (k0 + 2 * k1 + 2 * k2 + k3) / 6;             // Combined step
    vector yh = y + k * h;
    vector dy = (k - k0) * h;                               // error estimate

    return {yh, dy};
}
// Adaptive step size driver
std::pair<vector, std::vector<vector>> driver(
    const std::function<vector(double, const vector&)>& F, // Function dy/dx = f(x, y)
    std::pair<double, double> interval,                    // Interval [a, b]
    vector ystart,                                         // Initial conditions
    double h,                                              // Step size
    double acc,                                            // Desired accuracy
    double eps,                                            // Tolerance
    int nmax                                               // Maximum number of steps
) {
    double a = interval.first, b = interval.second;        
    double x = a;
    vector y = ystart.copy();
    vector xlist;
    xlist.push_back(x);
    std::vector<vector> ylist;
    ylist.push_back(y);
    
    int nsteps = 0; // Initialize step counter
    while(x < b && nsteps < nmax) {
        if (x + h > b) h = b - x;

        std::pair<vector, vector> step = rkstep4(F, x, y, h);
        vector yh = step.first;
        vector deltaY = step.second;

        double tol = (acc + eps * yh.norm()) * std::sqrt(h / (b - a));
        double err = deltaY.norm();
        
        if (err <= tol) {
            x += h;
            y = yh;
            xlist.push_back(x);
            ylist.push_back(y);
        }

        h = std::min(std::pow(tol / err, 0.25) * 0.95, 2.0) * h;
        nsteps++;
    }

    return {xlist, ylist}; 
}


// New driver function that includes interpolation
std::function<vector(double)> driver_interp(
    const std::function<vector(double, const vector&)>& F,
    std::pair<double, double> interval,
    vector ystart,
    double h,
    double acc,
    double eps,
    int nmax) {
    
    // Run the original ODE solver to generate xlist and ylist
    std::pair<vector, std::vector<vector>> solution = driver(F, interval, ystart, h, acc, eps, nmax);
    vector xlist = solution.first;
    std::vector<vector> ylist = solution.second;

    size_t dimensions = ystart.size;
    std::vector<CSpline> splines;

    vector xlist_custom(xlist);
    for (size_t dim = 0; dim < dimensions; ++dim) {
        std::vector<double> ys(ylist.size());
        for (size_t i = 0; i < ylist.size(); ++i) {
            ys[i] = ylist[i][dim];
        }
        vector ys_custom(ys);
        splines.emplace_back(xlist_custom, ys_custom);
    }

    std::function<vector(double)> interpolant = [splines, dimensions](double x) -> vector {
        vector y(dimensions);
        for (size_t dim = 0; dim < dimensions; ++dim) {
            y[dim] = splines[dim].evaluate(x);
        }
        return y;
    };

    return interpolant;
}