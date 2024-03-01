#include "../matrix/vector.h"
#include "../matrix/matrix.h"
#include "QRGS.h"
#include <cmath>
#include <iostream>


QRGS::QRGS(const matrix& M) : A(M.copy()), Q(M.copy()), R(matrix(M.size2, M.size2)) {
    size_t m = M.size2; // Number of columns

    // Initialization of R and modification of Q for QR decomposition
    for (size_t i = 0; i < m; ++i) {
        R(i, i) = Q[i].norm();
        Q.setCol(i, Q[i] / R(i, i));
        for (size_t j = i + 1; j < m; ++j) {
            R(i, j) = Q[i].dot(Q[j]);
            Q.setCol(j, Q[j] - Q[i] * R(i, j));
        }
    }
}


vector QRGS::solve(const vector& r) const {
    // Compute y = Q^Tb
    vector y = Q.transpose() * r; // Assuming matrix-vector multiplication is defined

    // Now solve Rx = y using back substitution
    vector x(R.size2);
    for (std::ptrdiff_t i = R.size2 - 1; i >= 0; i--) {
        double sum = 0;
        for (size_t k = i + 1; k < R.size2; k++) {
            sum += R(i, k) * x[k];
        }
        x[i] = (y[i] - sum) / R(i, i);
    }
    return x;
}

double QRGS::det() const {
    double determinant = 1.0;
    for (size_t i = 0; i < R.size1; ++i) {
        determinant *= R(i, i);
    }
    return determinant;
}

matrix QRGS::inverse() const {
    size_t m = A.size2;
    matrix B(m, m);
    vector e(A.size1), sol;
    for (size_t i = 0; i < m; ++i) {
        e[i] = 1;
        sol = solve(e);
        B.setCol(i, sol);
        e[i] = 0;
    }
    return B;
}