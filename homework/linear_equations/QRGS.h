#ifndef QRGS_H
#define QRGS_H
#include "../matrix/vector.h"
#include "../matrix/matrix.h"

struct QRGS {
    matrix A, Q, R;

    // Constructor to perform QR decomposition
    QRGS(const matrix& M);

    // Method to solve linear equations
    vector solve(const vector& r) const;

    // Method to calculate the determinant
    double det() const;

    // Method to calculate the inverse
    matrix inverse() const;
    
};

#endif // QRGS_H
