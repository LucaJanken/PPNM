#include "min.h"

// Numerical Gradient
vector num_grad(std::function<double(const vector&)> f, vector& x) {
    double eps = 1e-7;
    double fx = f(x);
    vector grad(x.size);
    for (size_t i = 0; i < x.size; i++) {
        double dx = abs(x[i]) * eps;
        if (abs(x[i]) < sqrt(eps)) {
            dx = eps;
        }
        x[i] += dx;
        grad[i] = (f(x) - fx) / dx;
        x[i] -= dx;
    }
    return grad;
}

// Backtracking Line Search with B_inv reset
double line_search(std::function<double(const vector&)> f, const vector& x, const vector& d, const vector& grad, matrix& B_inv) {
    double lambda = 1;
    while (f(x + lambda * d) > f(x) + 1e-4 * lambda * grad.dot(d)) {
        lambda *= 0.5;
        if (lambda < 1.0 / 1024) {
            B_inv = matrix::identity(B_inv.size1);
            break;
        }
    }
    return lambda;
}

// Symmetric Broyden's Update
void update(matrix& B_inv, const vector& s, const vector& y) {
    vector u = s - B_inv * y;
    double sy = s.dot(y);
    if (std::abs(sy) > 1e-6) {
        double gamma = (u.dot(y)) / (2.0 * sy);
        vector a = (u - gamma * s) / sy;
        for (size_t i = 0; i < B_inv.size1; ++i) {
            for (size_t j = 0; j < B_inv.size2; ++j) {
                B_inv(i, j) += a[i] * s[j] + s[i] * a[j];
            }
        }
    }
}

// SR1 Update
void sr1_update(matrix& B_inv, const vector& s, const vector& y) {
    vector v = s - B_inv * y;
    double vy = v.dot(y);

    if (std::abs(vy) > 1e-6) {
        matrix vvT(B_inv.size1, B_inv.size2);
        for (size_t i = 0; i < vvT.size1; ++i) {
            for (size_t j = 0; j < vvT.size2; ++j) {
                vvT(i, j) = v[i] * v[j];
            }
        }

        for (size_t i = 0; i < B_inv.size1; ++i) {
            for (size_t j = 0; j < B_inv.size2; ++j) {
                B_inv(i, j) += vvT(i, j) / vy; 
            }
        }
    }
}

// Quasi-Newton Method with adjustable update method
vector quasi_newton(std::function<double(const vector&)> f, vector x, double tol, bool useSR1) {
    size_t n = x.size;
    matrix B_inv = matrix::identity(n);

    vector grad = num_grad(f, x);

    size_t iter = 0;
    while (grad.norm() > tol) {
        vector d = B_inv * grad * -1.0;

        double lambda = line_search(f, x, d, grad, B_inv); 
        vector s = d * lambda; 
        vector x_new = x + s;
        vector y = num_grad(f, x_new) - grad;

        if (useSR1) {
            sr1_update(B_inv, s, y);
        } else {
            update(B_inv, s, y);
        }

        x = x_new;
        grad = num_grad(f, x);
        ++iter;
    }
    std::cout << "Quasi-Newton iterations: " << iter << std::endl;
    return x;
}