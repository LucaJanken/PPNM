#ifndef HAVE_MIN_H
#define HAVE_MIN_H

#include <functional>
#include <iostream>
#include "../matrix/vector.h"
#include "../matrix/matrix.h"

// Numerical Gradient
vector num_grad(std::function<double(const vector&)> f, const vector& x);

// Line Search
double line_search(std::function<double(const vector&)> f, const vector& x, const vector& d, const vector& grad, matrix& B_inv);

// Symmetric Broyden's Update
void update(matrix& B_inv, const vector& s, const vector& y);

// SR1 Update
void sr1_update(matrix& B_inv, const vector& s, const vector& y);

// Quasi-Newton Method
vector quasi_newton(std::function<double(const vector&)> f, vector x, double tol = 1e-6, bool useSR1 = false);

#endif // HAVE_MIN_H