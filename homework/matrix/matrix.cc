#include "matrix.h"
#include "vector.h"
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <iostream>

matrix::matrix() { // default constructor
    size1 = 0;
    size2 = 0;
    data = nullptr;
}

matrix::matrix(size_t n, size_t m) { // parameterized constructor
    size1 = n;
    size2 = m;
    data = new double[size1 * size2];
    for (size_t i = 0; i < size1 * size2; i++) data[i] = 0;
}

matrix::matrix(const matrix& other) { // copy constructor
    size1 = other.size1;
    size2 = other.size2;
    data = new double[size1 * size2];
    for (size_t i = 0; i < size1 * size2; i++) data[i] = other.data[i];
}

matrix::matrix(matrix&& tmp) { // move constructor
    size1 = tmp.size1;
    size2 = tmp.size2;
    data = tmp.data;
    tmp.size1 = 0;
    tmp.size2 = 0;
    tmp.data = nullptr;
}

matrix::~matrix() { // destructor
    size1 = 0; 
    size2 = 0; 
    delete[] data;
}

matrix& matrix::operator=(const matrix& other) { // copy assignment
    if (this == &other) return *this;
    (*this).~matrix();
    size1 = other.size1;
    size2 = other.size2;
    data = new double[size1 * size2];
    for (size_t i = 0; i < size1 * size2; i++) data[i] = other.data[i];
    return *this;
}

matrix& matrix::operator=(matrix&& tmp) { // move assignment
    (*this).~matrix();
    size1 = tmp.size1;
    size2 = tmp.size2;
    data = tmp.data;
    tmp.size1 = 0;
    tmp.size2 = 0;
    tmp.data = nullptr;
    return *this;
}

matrix& matrix::operator+=(const matrix& other) {
    assert(size1 == other.size1 && size2 == other.size2);
    for (size_t i = 0; i < size1 * size2; i++) data[i] += other.data[i];
    return *this;
}

matrix& matrix::operator-=(const matrix& other) {
    assert(size1 == other.size1 && size2 == other.size2);
    for (size_t i = 0; i < size1 * size2; i++) data[i] -= other.data[i];
    return *this;
}

matrix& matrix::operator*=(const matrix& other) {
    assert(size2 == other.size1);
    matrix R = matrix(size1, other.size2);
    for (size_t i=0; i<R.size1; i++) {
        for (size_t j=0; j<R.size2; j++) {
            double sum = 0;
            for (size_t k=0; k<size2; k++) sum += data[i * size2 + k] * other.data[k * other.size2 + j];
            R.data[i * R.size2 + j] = sum;
        }
    }
    (*this) = R;
    return *this;
}

matrix& matrix::operator*=(double x) {
    for (size_t i = 0; i < size1 * size2; i++) data[i] *= x;
    return *this;
}

matrix& matrix::operator/=(double x) {
    for (size_t i = 0; i < size1 * size2; i++) data[i] /= x;
    return *this;
}

matrix operator+(const matrix& A, const matrix& B) {
    assert(A.size1 == B.size1 && A.size2 == B.size2);
    matrix R = matrix(A.size1, A.size2);
    for (size_t i = 0; i < R.size1 * R.size2; i++) R.data[i] = A.data[i] + B.data[i];
    return R;
}

matrix operator-(const matrix& A, const matrix& B) {
    assert(A.size1 == B.size1 && A.size2 == B.size2);
    matrix R = matrix(A.size1, A.size2);
    for (size_t i = 0; i < R.size1 * R.size2; i++) R.data[i] = A.data[i] - B.data[i];
    return R;
}

matrix operator*(const matrix& A, const matrix& B) {
    assert(A.size2 == B.size1);
    matrix R = matrix(A.size1, B.size2);
    for (size_t k=0; k<A.size2; k++) {
        for (size_t j=0; j<B.size2; j++) {
            for (size_t i=0; i<A.size1; i++) {
                R.set(i, j, R(i, j) + A(i, k) * B(k, j));
            }
        }
    }
    return R;
}

matrix operator*(const matrix& A, double x) {
    matrix R = matrix(A.size1, A.size2);
    for (size_t i=0; i<R.size1 * R.size2; i++) R.data[i] = A.data[i] * x;
    return R;
}

matrix operator*(double x, const matrix& A) {
    matrix R = matrix(A.size1, A.size2);
    for (size_t i=0; i<R.size1 * R.size2; i++) R.data[i] = A.data[i] * x;
    return R;
}

matrix operator/(const matrix& A, double x) {
    matrix R = matrix(A.size1, A.size2);
    for (size_t i=0; i<R.size1 * R.size2; i++) R.data[i] = A.data[i] / x;
    return R;
}

vector operator*(const matrix& M, const vector& v) {
    assert(M.size2 == v.size);
    vector r(M.size1);
    for (size_t i=0; i<r.size; i++) {
        double sum = 0;
        for (size_t j=0; j<v.size; j++) sum += M(i, j) * v[j];
        r[i] = sum;
    }
    return r;   
}

void matrix::print(std::string s) const {
    printf("%s\n", s.c_str());
    for (size_t i=0; i<size1; i++) {
        for (size_t j=0; j<size2; j++) printf("%9.4g ", (*this)(i, j));
        printf("\n");
    }
}

// Generate an identity matrix of size n
matrix matrix::identity(size_t n) {
    matrix I(n, n);
    for (size_t i = 0; i < n; ++i) {
        I.data[i * n + i] = 1.0;
    }
    return I;
}

// Transpose the matrix
matrix matrix::transpose() const {
    matrix T(size2, size1);
    for (size_t i = 0; i < size1; ++i) {
        for (size_t j = 0; j < size2; ++j) {
            T.data[j + i * size2] = data[i + j * size1];
        }
    }
    return T;
}


// Compare two matrices with a tolerance
bool matrix::compare(const matrix& A, const matrix& B, double tol) {
    if (A.size1 != B.size1 || A.size2 != B.size2) return false;
    for (size_t i = 0; i < A.size1; ++i) {
        for (size_t j = 0; j < A.size2; ++j) {
            if (std::fabs(A(i, j) - B(i, j)) > tol) return false;
        }
    }
    return true;
}

// Check if the matrix is upper triangular
bool matrix::isUpTri(double tol) const {
    for (size_t i = 1; i < size1; ++i) {
        for (size_t j = 0; j < i && j < size2; ++j) {
            if (std::fabs(data[i + j * size1]) > tol) return false;
        }
    }
    return true;
}

// Check if the matrix is lower triangular
bool matrix::isLowTri(double tol) const {
    for (size_t i = 0; i < size1 - 1; ++i) {
        for (size_t j = i + 1; j < size2; ++j) {
            if (std::fabs(data[i + j * size1]) > tol) return false;
        }
    }
    return true;
}

// Set the i-th column of the matrix to a vector
void matrix::setCol(size_t i, const vector& v) {
    assert(v.size == size1);
    for (size_t j = 0; j < size1; ++j) {
        this->data[i * size1 + j] = v[j];
    }
}