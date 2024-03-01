#include "QRGS.h"
#include <iostream>
#include <random>

// Function to generate random matrix
matrix rnd_matrix(size_t n, size_t m) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    matrix A(n, m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            A.set(i, j, dis(gen));
        }
    }
    return A;
}

// Function to generate a random vector
vector rnd_vec(size_t n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    vector b(n);
    for (size_t i = 0; i <n; i++) {
        b[i] = dis(gen);
    }
    return b;
}

// Main function
int main() {

    // Task A
    std::cout << "QR decompisition: \n";
  
    // Generate random tall matrix
    size_t n = 8;
    size_t m = 5;
    matrix A = rnd_matrix(n, m);
    A.print("A:");

    // Perform QR factorisation
    QRGS qrgs(A);
    matrix Q = qrgs.Q;
    matrix R = qrgs.R;
    Q.print("Q:");
    R.print("R:");

    // Check for upper triangularity
    if (R.isUpTri()) {
        std::cout << "R is upper triangular. \n";
    } else {
        std::cout << "R is not upper triangular. \n";
    }

    // Ensure tat Q^T * Q = I
    matrix Q_T = Q.transpose();
    Q_T.print("Q^T:");
    matrix Q_TQ = Q_T * Q;
    Q_TQ.print("Q^T * Q:");
    matrix I1 = matrix::identity(Q.size2);
    if (matrix::compare(Q_TQ, I1)) {
        std::cout << "Q^T * Q is the identity matrix. \n";
    } else {
        std::cout << "Q^T * Q is not the identity matrix. \n";
    }

    // Ensure that Q * R = A
    matrix QR = Q * R;
    QR.print("Q * R:");
    if (matrix::compare(QR, A)) {
        std::cout << "Q * R equals A. \n";
    } else {
        std::cout << "Q * R does not equal A. \n";
    }

    std::cout << "QR solve: \n";

    // Generate square matrix and vector with random entries
    matrix A_square = rnd_matrix(m, m);
    A_square.print("A_square:");
    vector b = rnd_vec(m);
    b.print("b:");

    // Perform QR factorisation and solve
    QRGS qrgs_square(A_square);
    vector x = qrgs_square.solve(b);
    x.print("x:");

    // Ensure that A * x = b
    vector Ax = A_square * x;
    Ax.print("A * x:");
    if (compare(Ax, b)) {
        std::cout << "A * x equals b. \n";
    } else {
        std::cout << "A * x does not equal b. \n";
    }

    // Find determinant of A and A_square through their R's
    double det_A = qrgs.det();
    std::cout << "abs(det) of A via R: " << det_A << std::endl;

    double det_A_square = qrgs_square.det();
    std::cout << "abs(det) of A_square via R: " << det_A_square << std::endl;

    // Task B
    std::cout << "Matrix inverse";

    // Calculate Inverse
    matrix B = qrgs_square.inverse();
    B.print("B:");
    matrix AB = A_square * B;
    AB.print("A * B:");
    matrix I2 = matrix::identity(m);

    // Ensure A * B = I
    if (matrix::compare(AB, I2)) {
        std::cout << "AB equals I. \n";
    } else {
        std::cout << "AB does not equal I. \n";
    }

    return 0;
}
