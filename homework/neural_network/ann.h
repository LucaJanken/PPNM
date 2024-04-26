#ifndef ANN_H
#define ANN_H

#include "../matrix/vector.h"
#include "../minimization/min.h"
#include <functional>

class ANN {
private:
    int n;
    vector a, b, w;
    std::function<double(double)> activation;
    std::function<double(double)> activation_derivative;
    std::function<double(double)> activation_second_derivative;
    std::function<double(double)> activation_integral;

public:
    ANN(int n);
    double response(double x);
    double derivative(double x);
    double second_derivative(double x);
    double antiderivative(double x);  // Requires bounds or constant of integration if used absolutely
    void train(const vector& x, const vector& y);
};


#endif // ANN_H
