#ifndef HAVE_MATRIX_H
#define HAVE_MATRIX_H
#include <cstdlib>
#include <string>
#include "vector.h"

struct matrix {

    size_t size1, size2; double* data;

    matrix(size_t, size_t);         // parameterized constructor
    matrix(const matrix&);          // copy constructor
    matrix(matrix&&);               // move constructor
    matrix();                       // default constructor
    ~matrix();                      // destructor

    matrix& operator=(const matrix&);   // copy assignment
    matrix& operator=(matrix&&);        // move assignment

    matrix copy() const {
        matrix R(*this);
        return R;
    }

    void set (size_t i, size_t j, double x) {
        data[i + j * size1] = x;
    }

    double get (size_t i, size_t j) const {
        return data[i + j * size1];
    }

    double& operator()(size_t i, size_t j) {
        return data[i + j * size1];
    }

    double operator()(size_t i, size_t j) const {
        return data[i + j * size1];
    }

    vector operator[](size_t i) const {
        vector r(size1);
        for (size_t j = 0; j < size1; ++j) {
            r[j] = data[j + i * size1];
        }
        return r;
    }

    matrix& operator+=(const matrix&);
    matrix& operator-=(const matrix&);
    matrix& operator*=(const matrix&);
    matrix& operator*=(double);
    matrix& operator/=(double);
    matrix& operator^(int);

    static matrix identity(size_t n); // Generate identity matrix
    matrix transpose() const;
    static bool compare(const matrix& A, const matrix& B, double tol = 1e-9);
    bool isUpTri(double tol = 1e-9) const;
    bool isLowTri(double tol = 1e-9) const;
    void setCol(size_t, const vector&);
    vector getCol(size_t) const;

    void print(std::string s="") const;

};

matrix operator+(const matrix&, const matrix&);
matrix operator-(const matrix&, const matrix&);
matrix operator*(const matrix&, const matrix&);
matrix operator*(const matrix&, double x);
matrix operator*(double x, const matrix&);
matrix operator/(const matrix&, double x);
vector operator*(const matrix&, const vector&);


#endif