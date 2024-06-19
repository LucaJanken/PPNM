#ifndef RAYLEIGH_QUOTIENT_H
#define RAYLEIGH_QUOTIENT_H

#include "../homework/matrix/vector.h"
#include "../homework/matrix/matrix.h"

double RayleighQuotient(const vector& v, const matrix& H);
vector Gradient(const vector& v, const matrix& H, double R);

#endif // RAYLEIGH_QUOTIENT_H
