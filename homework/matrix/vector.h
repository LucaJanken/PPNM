#ifndef HAVE_VECTOR_H
#define HAVE_VECTOR_H
#include <cstdlib>
#include <cassert>

struct vector {
    
    size_t size; double* data;

    vector(size_t);         // parameterized constructor
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

    double dot(const vector&) const;
    double norm() const;


    vector& operator+=(const vector&);
    vector& operator-=(const vector&);
    vector& operator*=(const vector&);
    vector& operator*=(double);
    vector& operator/=(double);

    void print(const char* s="") const;

};

vector operator+(const vector&, const vector&);
vector operator-(const vector&, const vector&);
vector operator*(const vector&, const vector&);
vector operator*(const vector&, double x);
vector operator*(double x, const vector&);
vector operator/(const vector&, double x);

bool compare(const vector& a, const vector& b, double tol = 1e-9);

#endif