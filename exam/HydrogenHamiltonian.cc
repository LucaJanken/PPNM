#include "HydrogenHamiltonian.h"

matrix Hamiltonian(int rmax, double dr) {
    int npoints = static_cast<int>(rmax / dr) - 1;
    vector r(npoints);
    for (int i = 0; i < npoints; i++) {
        r[i] = dr * (i + 1);
    }
    matrix H(npoints, npoints);
    for (int i = 0; i < npoints - 1; i++) {
        H.set(i, i, -2);
        H.set(i, i + 1, 1);
        H.set(i + 1, i, 1); 
    }
    H.set(npoints - 1, npoints - 1, -2);
    H = H * (-0.5 / dr / dr);
    for (int i = 0; i < npoints; i++) {
        H.set(i, i, H.get(i, i) - 1 / r[i]);
    }
    return H;
}
