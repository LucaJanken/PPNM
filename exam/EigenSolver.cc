#include "EigenSolver.h"
#include "RayleighQuotient.h"
#include <cmath>
#include <limits>
#include <iostream>
#include <chrono>

void Normalize(vector& v) {
    double norm = 0.0;
    for (double vi : v) {
        norm += vi * vi;
    }
    norm = std::sqrt(norm);
    for (double& vi : v) {
        vi /= norm;
    }
}

double line_search(const std::function<double(const vector&)>& f, const vector& v, const vector& delta_v) {
    double alpha = 1.0;
    double rho = 0.5;
    double c = 1e-4;
    while (f(v - alpha * delta_v) > f(v) - c * alpha * delta_v.dot(delta_v)) {
        alpha *= rho;
    }
    return alpha;
}

vector newton_method(const std::function<double(const vector&)>& f, const std::function<vector(const vector&)>& grad_f, vector v, double tol, const std::chrono::time_point<std::chrono::high_resolution_clock>& start_time, double timeout) {
    size_t iter = 0;

    while (true) {
        auto now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = now - start_time;
        if (elapsed.count() > timeout) {
            throw std::runtime_error("Timeout exceeded");
        }

        vector grad = grad_f(v);
        if (grad.norm() < tol) break;

        vector delta_v = grad;
        double alpha = line_search(f, v, delta_v);
        v = v - alpha * delta_v;
        Normalize(v);

        ++iter;
    }

    std::cout << "Newton's method iterations: " << iter << std::endl;
    return v;
}

std::pair<double, vector> FindEigenPair(const matrix& H, double tol, double timeout) {
    int n = H.size1;
    vector v(n, 1.0);
    Normalize(v);

    auto start_time = std::chrono::high_resolution_clock::now();

    auto RayleighQuotientFunc = [&H](const vector& x) {
        return RayleighQuotient(x, H);
    };

    auto GradientFunc = [&H](const vector& x) {
        double R = RayleighQuotient(x, H);
        return Gradient(x, H, R);
    };

    v = newton_method(RayleighQuotientFunc, GradientFunc, v, tol, start_time, timeout);
    double eigenvalue = RayleighQuotient(v, H);

    return {eigenvalue, v};
}

std::pair<double, vector> FindSecondEigenPair(const matrix& H, const vector& first_eigenvector, double tol, double timeout) {
    int n = H.size1;

    // Convert the eigenvector into a column matrix
    matrix v1_column(n, 1);
    for (int i = 0; i < n; ++i) {
        v1_column.set(i, 0, first_eigenvector[i]);
    }

    // Convert the eigenvector into a row matrix
    matrix v1_row = v1_column.transpose();

    // Compute the outer product: v1_column * v1_row
    matrix outer_product = v1_column * v1_row;

    // Construct the projection matrix: P = I - outer_product
    matrix P = matrix::identity(n) - outer_product;

    // Project the original matrix: H' = P * H * P
    matrix H_prime = P * H * P;

    // Find the smallest eigenvalue of the modified matrix H_prime
    return FindEigenPair(H_prime, tol, timeout);
}
