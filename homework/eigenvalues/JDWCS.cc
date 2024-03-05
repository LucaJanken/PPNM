#include <cmath>
#include "JDWCS.h"
#include "../matrix/vector.h"
#include "../matrix/matrix.h"

void JDWCS::timesJ(matrix& A, int p, int q, double theta) {
    double c = cos(theta), s = sin(theta);
    for (size_t i = 0; i < A.size1; ++i) {
        double aip = A(i, p), aiq = A(i, q);
        A.set(i, p, c * aip - s * aiq);
        A.set(i, q, s * aip + c * aiq);
    }
}

void JDWCS::Jtimes(matrix& A, int p, int q, double theta) {
    double c = cos(theta), s = sin(theta);
    for (size_t j = 0; j < A.size1; ++j) {
        double apj = A(p, j), aqj = A(q, j);
        A.set(p, j, c * apj - s * aqj);
        A.set(q, j, s * apj + c * aqj);
    }
}

std::pair<matrix, matrix> JDWCS::cyclic(const matrix& A) {
    int n = A.size1;
    matrix D = A.copy();
    matrix V = matrix::identity(n); // Initialize V as an identity matrix of size n
    int sweeps = 0;
    bool changed;

    do {
        sweeps++;
        changed = false;
        for (int q = n - 1; q > 0; q--) {
            for (int p = 0; p < q; p++) {
                double dpp = D.get(p, p);
                double dqq = D.get(q, q);
                double dpq = D.get(p, q);
                double phi = 0.5 * std::atan2(2 * dpq, dqq - dpp);
                double c = std::cos(phi), s = std::sin(phi);

                double dpp1 = c * c * dpp - 2 * s * c * dpq + s * s * dqq;
                double dqq1 = s * s * dpp + 2 * s * c * dpq + c * c * dqq;
                if (dpp1 != dpp || dqq1 != dqq) {
                    changed = true;
                    D.set(p, p, dpp1);
                    D.set(q, q, dqq1);
                    D.set(p, q, 0.0);
                    D.set(q, p, 0.0);

                    // Additional rotations for elements of D
                    for (int i = 0; i < p; i++) {
                        double dip = D.get(i, p), diq = D.get(i, q);
                        D.set(i, p, c * dip - s * diq);
                        D.set(i, q, s * dip + c * diq);
                        D.set(p, i, D.get(i, p));
                        D.set(q, i, D.get(i, q));
                    }
                    for (int i = p + 1; i < q; i++) {
                        double dpi = D.get(p, i), diq = D.get(i, q);
                        D.set(p, i, c * dpi - s * diq);
                        D.set(i, q, s * dpi + c * diq);
                        D.set(i, p, D.get(p, i));
                        D.set(q, i, D.get(i, q));
                    }
                    for (int i = q + 1; i < n; i++) {
                        double dpi = D.get(p, i), dqi = D.get(q, i);
                        D.set(p, i, c * dpi - s * dqi);
                        D.set(q, i, s * dpi + c * dqi);
                        D.set(i, p, D.get(p, i));
                        D.set(i, q, D.get(q, i));
                    }

                    // Update V to accumulate transformations
                    for (int i = 0; i < n; i++) {
                        double vip = V.get(i, p), viq = V.get(i, q);
                        V.set(i, p, c * vip - s * viq);
                        V.set(i, q, s * vip + c * viq);
                    }
                }
            }
        }
    } while (changed);

    return std::pair<matrix, matrix>(D, V);
}
