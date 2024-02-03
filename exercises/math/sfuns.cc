#include <cmath>
#include <limits> // Include for std::numeric_limits

double my_gamma(double x) {
    //single precision gamma function (formula from Wikipedia)
    if (x < 0) return M_PI / sin(M_PI * x) / my_gamma(1 - x); // Euler's reflection formula
    if (x < 9) return my_gamma(x + 1) / x; // Recurrence relation
    double lngamma = x * log(x + 1 / (12 * x - 1 / x / 10)) - x + log(2 * M_PI / x) / 2;
    return exp(lngamma);
}

double my_lngamma(double x) {
    if (x <= 0) return std::numeric_limits<double>::quiet_NaN();
    if (x == 1 || x == 2) return 0; // Avoid floating-point imprecision for ln(1)
    if (x < 9) return my_lngamma(x + 1) - log(x);
    double lngamma = x * log(x + 1 / (12 * x - 1 / x / 10)) - x + log(2 * M_PI / x) / 2;
    return lngamma;
}
