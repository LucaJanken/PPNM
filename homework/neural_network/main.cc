#include "ann.h"
#include <fstream>
#include <iostream>
#include <cmath>

// Define the derivatives and antiderivative of the target function
double target_function(double x) {
    return std::cos(5 * x - 1) * std::exp(-x * x);
}

int main() {
    const int num_points = 30;
    vector x_values(num_points), y_values(num_points);

    // Notify the sampling process and the number of points
    std::cout << "Sampling target function over " << num_points 
              << " points and training the neural network to generate responses.\n";

    for (size_t i = 0; i < num_points; i++) {
        double interval = 2.0 / (num_points - 1);
        double x = -1 + i * interval; // interval from -1 to 1
        x_values.set(i, x);
        y_values.set(i, target_function(x));
    }

    ANN neural_network(5); // using 5 neurons in the hidden layer
    neural_network.train(x_values, y_values);

    std::ofstream file("ann_results.csv");
    // file << "x, Target Function, ANN Response, Target First Derivative, ANN First Derivative, Target Second Derivative, ANN Second Derivative, ANN Antiderivative\n";
    for (size_t i = 0; i < x_values.size; i++) {
        double x = x_values.get(i);
        file << x << " "
             << target_function(x) << " "
             << neural_network.response(x) << " "
             << neural_network.derivative(x) << " "
             << neural_network.second_derivative(x) << " "
             << neural_network.antiderivative(x) << "\n";

    }
    file.close();

    std::cout << "ANN training complete. Results saved to ann_results.csv." << std::endl;
    return 0;
}
