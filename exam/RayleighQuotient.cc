#include "RayleighQuotient.h"
#include <cmath>

double RayleighQuotient(const vector& v, const matrix& H) {
    double numerator = 0.0, denominator = 0.0;
    int n = v.size;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            numerator += v[i] * H[i][j] * v[j];
        }
        denominator += v[i] * v[i];
    }
    return numerator / denominator;
}

vector Gradient(const vector& v, const matrix& H, double R) {
    int n = v.size;
    vector grad(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            grad[i] += H[i][j] * v[j];
        }
        grad[i] -= R * v[i];
    }
    return grad;
}
