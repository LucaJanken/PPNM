#include "vector.h"
#include <algorithm>
#include <cstdlib>
#include <cassert>
#include <cmath>

vector::vector() { // default constructor
    size = 0;
    data = nullptr;
}

vector::vector(size_t n, double init_value /* = 0.0 */): size(n), data(new double[n]) {
    for (size_t i = 0; i < size; ++i) {
        data[i] = init_value;
    }
}

vector::vector(const vector& other) { // copy constructor
    size = other.size;
    data = new double[size];
    for (size_t i=0; i<(*this).size; i++) data[i] = other.data[i];
}

vector::vector(vector&& tmp) { // move constructor
    size = tmp.size;
    data = tmp.data;
    tmp.size = 0;
    tmp.data = nullptr;
}

vector::~vector() { // destructor
    size = 0; delete[] data;
}

vector& vector::operator=(const vector& other) { // copy assignment
    if (this == &other) return *this;
    (*this).~vector();
    size = other.size;
    data = new double[size];
    for (size_t i=0; i<(*this).size; i++) data[i] = other.data[i];
    return *this;
}

vector& vector::operator=(vector&& tmp) { // move assignment
    (*this).~vector();
    size = tmp.size;
    data = tmp.data;
    tmp.size = 0;
    tmp.data = nullptr;
    return *this;
}

vector& vector::operator+=(const vector& other) { // addition assignment
    assert(size == other.size);
    for (size_t i=0; i<(*this).size; i++) data[i] += other.data[i];
    return *this;
}

vector& vector::operator-=(const vector& other) { // subtraction assignment
    assert(size == other.size);
    for (size_t i=0; i<(*this).size; i++) data[i] -= other.data[i];
    return *this;
}

vector& vector::operator*=(const vector& other) { // dot product assignment
    assert(size == other.size);
    for (size_t i=0; i<(*this).size; i++) data[i] *= other.data[i];
    return *this;
}

vector& vector::operator*=(double x) { // scalar multiplication assignment
    for (size_t i=0; i<(*this).size; i++) data[i] *= x;
    return *this;
}

vector& vector::operator/=(double x) { // scalar division assignment
    for (size_t i=0; i<(*this).size; i++) data[i] /= x;
    return *this;
}

vector operator/(const vector& v, double x) {
    vector result = vector(v.size);
    for (size_t i=0; i<result.size; i++) result.data[i] = v.data[i] / x;
    return result;
}

vector operator*(const vector& a, const vector& b) {
    assert(a.size == b.size);
    vector result = vector(a.size);
    for (size_t i=0; i<result.size; i++) result.data[i] = a.data[i] * b.data[i];
    return result;
}

vector operator*(const vector& a, double x) {
    vector result = vector(a.size);
    for (size_t i=0; i<result.size; i++) result.data[i] = a.data[i] * x;
    return result;
}

vector operator*(double x, const vector& a) {
    vector result = vector(a.size);
    for (size_t i=0; i<result.size; i++) result.data[i] = a.data[i] * x;
    return result;
}

vector operator+(const vector& a, const vector& b) {
    assert(a.size == b.size);
    vector result = vector(a.size);
    for (size_t i=0; i<result.size; i++) result.data[i] = a.data[i] + b.data[i];
    return result;
}

vector operator-(const vector& a, const vector& b) {
    assert(a.size == b.size);
    vector result = vector(a.size);
    for (size_t i=0; i<result.size; i++) result.data[i] = a.data[i] - b.data[i];
    return result;
}

void vector::print(const char* s) const {
    printf("%s\n", s);
    for (size_t i=0; i<size; i++) printf("%9.4g ", (*this)[i]);
    printf("\n");
}

void vector::resize(size_t new_size, double init_value) {
    double* new_data = new double[new_size];
    for (size_t i = 0; i < std::min(size, new_size); ++i) {
        new_data[i] = data[i];
    }
    for (size_t i = size; i < new_size; ++i) {
        new_data[i] = init_value;
    }
    delete[] data;
    data = new_data;
    size = new_size;
}

bool compare(const vector& a, const vector& b, double tol) {
    if (a.size != b.size) return false;

    for (size_t i = 0; i < a.size; ++i) {
        if (std::fabs(a.data[i] - b.data[i]) > tol) return false;
    }

    return true;
}

double vector::dot(const vector& other) const {
    assert(size == other.size);
    double sum = 0;
    for (size_t i=0; i<size; i++) sum += data[i] * other.data[i];
    return sum;
}

double vector::norm() const {
    double sum = 0;
    for (size_t i=0; i<size; i++) sum += data[i] * data[i];
    return sqrt(sum);
}
