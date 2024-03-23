#ifndef HAVE_MONTECARLO_H
#define HAVE_MONTECARLO_H

#include <cmath>
#include <random>
#include <functional>
#include <iostream>
#include "../matrix/vector.h"
#include "../matrix/matrix.h"

std::pair<double, double> plainmc(std::function<double(vector&)> f, const vector& a, const vector& b, int N);

std::pair<double, double> quasi_mc(std::function<double(vector&)> f, const vector& a, const vector& b, int N);

std::pair<double, double> stratified_mc(std::function<double(vector&)> f, const vector& a, const vector& b, int N, int nmin);

#endif