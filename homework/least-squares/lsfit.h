#ifndef HAVE_LS_FIT_H
#define HAVE_LS_FIT_H
#include "../matrix/vector.h"
#include "../matrix/matrix.h"
#include <functional>
#include <vector>
#include <tuple>

std::tuple<vector, matrix> lsfit(const std::vector<std::function<double(double)>>& fs, const vector& x, const vector& y, const vector& dy);

#endif // HAVE_LS_FIT_H