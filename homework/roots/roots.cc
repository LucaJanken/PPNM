#include "roots.h"

matrix calcJac(vfunc f, vector x) {
    vector f0 = f(x);
    int dim = x.size, outDim = f0.size;
    matrix J = matrix(outDim, dim);
    double eps = pow(2, -26);
    vector xPerturb = x.copy(), df(outDim);
    
    for (int i = 0; i < dim; i++) {
        double dx = eps * (abs(xPerturb[i]) + 1e-4);
        xPerturb[i] += dx;
        df = f(xPerturb) - f0;
        J.setCol(i, df / dx);
        xPerturb[i] = x[i];
    }
    return J;
}

vector newton(vfunc f, vector x, double eps) {
    double lambda = 1;
    vector x0 = x.copy(), Dx = x.copy();
    
    while (f(x0).norm() >= eps && Dx.norm() >= pow(2, -26) * x0.norm()) {
        vector f0 = f(x0);
        matrix J = calcJac(f, x0);
        QRGS solver(J);
        Dx = solver.solve(-f0);
        lambda = 1;
        
        while (f(x0 + lambda * Dx).norm() > (1 - lambda / 2) * f0.norm() && lambda > 1 / 1024) {
            lambda /= 2;
        }
        x0 += lambda * Dx;
    }
    return x0;
}