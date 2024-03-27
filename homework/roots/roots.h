#ifndef HAVE_ROOTS_H
#define HAVE_ROOTS_H

#include "../matrix/vector.h"
#include "../matrix/matrix.h"
#include "../linear_equations/QRGS.h"

using vfunc = std::function<vector(const vector&)>;

vector newton(vfunc f, vector x, double eps = 1e-3);

#endif

