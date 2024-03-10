#include <cmath>
#include <fstream>
#include "../matrix/vector.h"
#include "spline.h"

int main() {

    // Linear spline test
    int n = 10;
    vector x(n);
    vector y(n);

    for (size_t i = 0; i < x.size; i++) {
        x[i] = i;
        y[i] = std::cos(x[i]);
    }

    LSpline ls(x, y);

    std::ofstream file("lspline.csv");
    file << "x, f(x), F(x)\n";

    for (double xi = 0; xi < n - 1; xi += 0.1) {
        file << xi << ", " << ls.evaluate(xi) << ", " << ls.integral(xi) << "\n";
    }

    file.close();

    // Quadratic spline test
    vector x2(n);
    vector y2(n);

    for (size_t i = 0; i < x2.size; i++) {
        x2[i] = i;
        y2[i] = std::sin(x2[i]);
    }

    QSpline qs(x2, y2);

    std::ofstream file2("qspline.csv");
    file2 << "x, f(x), F(x), df/dx\n";

    for (double xi = 0; xi < n - 1; xi += 0.1) {
        file2 << xi << ", " << qs.evaluate(xi) << ", " << qs.integral(xi) << ", " << qs.derivative(xi) << "\n";
    }

    file2.close();
    
    // Cubic spline test
    vector x3(n);
    vector y3(n);

    for (size_t i = 0; i < x3.size; i++) {
        x3[i] = i;
        y3[i] = std::sin(x3[i]);
    }

    CSpline cs(x3, y3);
    
    std::ofstream file3("cspline.csv");
    file3 << "x, f(x), F(x), df/dx\n";

    for (double xi = 0; xi < n - 1; xi += 0.1) {
        file3 << xi << ", " << cs.evaluate(xi) << ", " << cs.integral(xi) << ", " << cs.derivative(xi) << "\n";
    }

    file3.close();
    
    return 0;
}