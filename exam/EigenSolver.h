#ifndef EIGEN_SOLVER_H
#define EIGEN_SOLVER_H

#include "../homework/matrix/vector.h"
#include "../homework/matrix/matrix.h"
#include <utility>
#include <chrono>
#include <functional>

void Normalize(vector& v);
double line_search(const std::function<double(const vector&)>& f, const vector& v, const vector& delta_v);
vector newton_method(const std::function<double(const vector&)>& f, const std::function<vector(const vector&)>& grad_f, vector v, double tol, const std::chrono::time_point<std::chrono::high_resolution_clock>& start_time, double timeout);
std::pair<double, vector> FindEigenPair(const matrix& H, double tol = 1e-10, double timeout = 5.0);
std::pair<double, vector> FindSecondEigenPair(const matrix& H, const vector& first_eigenvector, double tol = 1e-10, double timeout = 5.0);

#endif // EIGEN_SOLVER_H
