#include <iostream>
#include <limits>
#include <cmath>
#include "math_utils.h"


int find_max_int() {
    int i = 1;
    while (i + 1 > 1) {
        i = i + 1;
    }
    return i;
}


int find_min_int() {
    int i = 1;
    while (i - 1 < 1) {
        i = i - 1;
    }
    return i;
}


double double_epsilon() {
    double x = 1.0;
    while (1.0 + x != 1.0) {
        x = x / 2.0;
    }
    return x * 2.0;
}


float float_epsilon() {
    float y = 1.0f;
    while ((float)(1.0f + y) != 1.0f) {
        y = y / 2.0f;
    }
    return y * 2.0f;
}


double sum_method_A(int n, double epsilon) {
    double tiny = epsilon / 2.0;
    double sumA = 1;
    for (int i = 0; i < n; i = i + 1) {
        sumA = sumA + tiny;
    }
    return sumA;
}

double sum_method_B(int n, double epsilon) {
    double tiny = epsilon / 2.0;
    double sumB = 0;
    for (int i = 0; i < n; i = i + 1) {
        sumB = sumB + tiny;
    }
    return sumB + 1;
}

bool approx(double a, double b, double acc, double eps) {
    if (std::abs(a - b) <= acc) return true;
    if (std::abs(a - b) <= std::max(std::abs(a), std::abs(b)) * eps) return true;
    return false;
}