#ifndef HAVE_SPLINE_H
#define HAVE_SPLINE_H
#include "../matrix/vector.h"

// Abstract base class for a generic spline
class Spline {
public:
    // Constructor
    Spline() {}

    // Virtual destructor to ensure proper cleanup of derived types
    virtual ~Spline() {}

    // Pure virtual methods for interpolation and integration to be implemented by derived classes
    virtual double evaluate(double z) const = 0;
    virtual double derivative(double z) const = 0;
    virtual double integral(double z) const = 0;

protected:
    // Utility method for binary search, useful for all spline types
    static int binsearch(const vector& x, double z);
};

// Derived class for linear spline interpolation
class LSpline : public Spline {
private:
    vector x, y;

public:
    LSpline(const vector& xs, const vector& ys);
    virtual ~LSpline() {}

    virtual double evaluate(double z) const override;
    virtual double derivative(double z) const override { (void)z; return 0; }
    virtual double integral(double z) const override;
};

// Derived class for quadratic spline interpolation
class QSpline : public Spline {
private:
    vector x, y, b, c;

public:
    QSpline(const vector& xs, const vector& ys);
    virtual ~QSpline() {}

    virtual double evaluate(double z) const override;
    virtual double derivative(double z) const override;
    virtual double integral(double z) const override;
};

// Derived class for cubic spline interpolation
class CSpline : public Spline {
private:
    vector x, y, b, c, d;

public:
    CSpline(const vector& xs, const vector& ys);
    virtual ~CSpline() {}

    virtual double evaluate(double z) const override;
    virtual double derivative(double z) const override;
    virtual double integral(double z) const override;
};

#endif // HAVE_SPLINE_H
