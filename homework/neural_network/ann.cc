#include "ann.h"
#include <cmath>

ANN::ANN(int num_neurons) : n(num_neurons), a(num_neurons), b(num_neurons), w(num_neurons) {
    double start = -1.0, end = 1.0, interval = (end - start) / (n - 1);
    for (int i = 0; i < n; i++) {
        a.set(i, start + i * interval);
        b.set(i, 0.5);
        w.set(i, 0.5);
    }
    activation = [](double x) { return x * std::exp(-x * x); };
    activation_derivative = [](double x) { return (1 - 2 * x * x) * std::exp(-x * x); };
    activation_second_derivative = [](double x) { return (4 * x * x * x - 6 * x) * std::exp(-x * x); };
    activation_integral = [](double x) { return -0.5 * std::exp(-x * x); };
}

double ANN::response(double x) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size; ++i) {
        double transformed_input = (x - a[i]) / b[i];
        sum += activation(transformed_input) * w[i];
    }
    return sum;
}

double ANN::derivative(double x) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size; ++i) {
        double transformed_input = (x - a[i]) / b[i];
        sum += activation_derivative(transformed_input) * w[i] / b[i];
    }
    return sum;
}

double ANN::second_derivative(double x) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size; ++i) {
        double transformed_input = (x - a[i]) / b[i];
        sum += activation_second_derivative(transformed_input) * w[i] / (b[i] * b[i]);
    }
    return sum;
}

double ANN::antiderivative(double x) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size; ++i) {
        double transformed_input = (x - a[i]) / b[i];
        sum += activation_integral(transformed_input) * w[i] * b[i];
    }
    return sum;
}

void ANN::train(const vector& x, const vector& y) {
    vector initial_params(3 * n);
    for (int i = 0; i < n; i++) {
        initial_params.set(i, a.get(i));
        initial_params.set(n + i, b.get(i));
        initial_params.set(2 * n + i, w.get(i));
    }

    auto cost_function = [this, &x, &y](const vector& params) {
        for (int i = 0; i < n; i++) {
            a.set(i, params.get(i));
            b.set(i, params.get(n + i));
            w.set(i, params.get(2 * n + i));
        }

        double cost = 0.0;
        for (size_t k = 0; k < x.size; ++k) {
            double error = response(x.get(k)) - y.get(k);
            cost += error * error;
        }
        return cost;
    };

    vector optimized_params = quasi_newton(cost_function, initial_params, 1e-3);

    for (int i = 0; i < n; i++) {
        a.set(i, optimized_params.get(i));
        b.set(i, optimized_params.get(n + i));
        w.set(i, optimized_params.get(2 * n + i));
    }
}
