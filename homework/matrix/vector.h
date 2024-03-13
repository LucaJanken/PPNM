#ifndef HAVE_VECTOR_H
#define HAVE_VECTOR_H
#include <cstdlib>
#include <cassert>
#include <vector>
#include <algorithm>

struct vector {
    
    size_t size; double* data;

    vector(size_t, double init_value = 0.0);         // parameterized constructor
    vector(const vector&);  // copy constructor
    vector(vector&&);       // move constructor
    vector();               // default constructor
    ~vector();              // destructor

    vector& operator=(const vector&);   // copy assignment
    vector& operator=(vector&&);        // move assignment

    vector copy() {
        vector r(*this); 
        return r;
    }

    double& operator[](size_t i) {
        assert(i < size);
        return data[i];
    }

    const double& operator[](size_t i) const {
        assert(i < size);
        return data[i];
    }

    void set(size_t i, double x) {
        assert(i < size);
        data[i] = x;
    }
    
    double get(size_t i) const {
        assert(i < size);
        return data[i];
    }

    double* begin() {
        return data;
    }

    double* end() {
        return data + size;
    }

    void push_back(double x) {
        resize(size + 1, x);
    }
         

    double dot(const vector&) const;
    double norm() const;


    vector& operator+=(const vector&);
    vector& operator-=(const vector&);
    vector& operator*=(const vector&);
    vector& operator*=(double);
    vector& operator/=(double);

    void print(const char* s="") const;
    void resize(size_t new_size, double init_value = 0.0);

    // Constructor that accepts a std::vector<double>
    vector(const std::vector<double>& init_values);

};

vector operator+(const vector&, const vector&);
vector operator-(const vector&, const vector&);
vector operator*(const vector&, const vector&);
vector operator*(const vector&, double x);
vector operator*(double x, const vector&);
vector operator/(const vector&, double x);

bool compare(const vector& a, const vector& b, double tol = 1e-9);

#endif