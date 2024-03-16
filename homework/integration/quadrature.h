#ifndef HAVE_QUADRATURE_H
#define HAVE_QUADRATURE_H
#include <functional>
#include <cmath>
#include <tuple>
#include <limits>

std::tuple<double, double, int> integrate(const std::function<double(double)>& f, double a, double b, double delta = 0.001, double eps = 0.001, double f2 = NAN, double f3 = NAN);

std::tuple<double, double, int> transint(const std::function<double(double)>& f, double a, double b, double delta, double eps);

std::tuple<double, double, int> integrate_infinite(const std::function<double(double)>& f, double a, double b, double delta = 0.001, double eps = 0.001);

#endif