#include "vector.h"
#include "matrix.h"
#include <iostream>

int main() {
    // Initialize matrix dimensions
    size_t n = 3, m = 3;

    // Create and populate matrix A
    matrix A(n, m);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < m; j++) {
            A(i, j) = i + j + 1; // Assuming this meant to use [] instead of , for access
        }
    }
    A.print("A=");

    // Copy matrix A to B and print
    matrix B = A;
    B.print("B=");

    // Modify and print matrix B
    B(1, 0) = 9; // Correcting index access syntax
    B.print("B after B(0, 0)=9");

    // Print A after B modification to verify immutability
    A.print("A after B(0, 0)=9");

    // Operations with A and printing results
    (A * 2).print("A*2=");
    (A / 2).print("A/2=");

    // Additional matrix operation
    matrix C = A * 2;
    C.print("C=");

    // Create and populate vector v
    vector v(m);
    for (size_t i = 0; i < v.size; i++) { // Assuming a size() method exists
        v[i] = 2 * i + 1;
    }
    v.print("v=");

    // Vector operations and printing results
    (v + v).print("v+v=");
    (v - v).print("v-v=");
    (v * 2).print("v*2=");

    // Set colum 2 of matrix A to vector v and print
    A.setCol(2, v);
    A.print("A after A.setCol(2, v)");

    // Test matrix matrix multiplication
    matrix D = A * B;
    matrix NI = -1 * matrix::identity(3);
    A.print("A=");
    B.print("B=");
    D.print("D=");

    matrix H(3, 3);
    H.set(0, 0, 0.71); H.set(0, 1, 0.4053); H.set(0, 2, 0.9491);
    H.set(1, 0, 0.8562); H.set(1, 1, 0.9863); H.set(1, 2, 0.134);
    H.set(2, 0, 0.9444); H.set(2, 1, 0.1619); H.set(2, 2, 0.3816);

    matrix H_inv(3, 3);
    H_inv.set(0, 0, -0.609623394015216); H_inv.set(0, 1, 0.00172429355436418); H_inv.set(0, 2, 1.51562502076404);
    H_inv.set(1, 0, 0.344065170418702); H_inv.set(1, 1, 1.07493387305472); H_inv.set(1, 2, -1.23321119558103);
    H_inv.set(2, 0, 1.36274680874524); H_inv.set(2, 1, -0.460325515933702); H_inv.set(2, 2, -0.607178661019358);

    //Identity test
    matrix I_test = H * H_inv;
    I_test.print("I_test=");
    return 0;
}

