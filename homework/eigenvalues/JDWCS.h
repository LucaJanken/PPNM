#ifndef HAVE_JDWCS_H
#define HAVE_JDWCS_H
#include "../matrix/vector.h"
#include "../matrix/matrix.h"
// Struct to perform Jacobi diagonalization with cyclic sweeps

struct JDWCS {

    static void timesJ(matrix& A, int p, int q, double theta);

    static void Jtimes(matrix& A, int p, int q, double theta);

    static std::pair<matrix, matrix> cyclic(const matrix& M);

};

#endif // HAVE_JDWCS_H