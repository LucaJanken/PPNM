#include "../matrix/vector.h"
#include "../matrix/matrix.h"
#include "../linear_equations/QRGS.h"
#include <functional>
#include <vector>
#include <tuple>

std::tuple<vector, matrix> lsfit(const std::vector<std::function<double(double)>>& fs, const vector& x, const vector& y, const vector& dy) {
    if (x.size != y.size || x.size != dy.size || y.size != dy.size) {
        throw std::invalid_argument("x, y, and dy must have the same size");
    }
        matrix A = matrix(x.size, fs.size());
        vector b = vector(y.size);
        for (size_t i = 0; i < x.size; i++) {
            b.set(i, y.get(i) / dy.get(i));
            for (size_t k = 0; k < fs.size(); k++) {
                A.set(i, k, fs[k](x.get(i)) / dy.get(i));
            }
        }
        QRGS qrgs(A);
        matrix Q = qrgs.Q;
        matrix R = qrgs.R;
        vector c = qrgs.solve(b);
        matrix Rinv = qrgs.inverse();
        matrix cov = Rinv * Rinv.transpose();
        return std::make_tuple(c, cov);
}
