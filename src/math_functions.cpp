/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "math_functions.hpp"

#include <cmath>

namespace netreg
{

    double abs_sum(matrix<double> &source1, matrix<double> &source2)
    {
        const int ncol = static_cast<int>(source1.n_cols);
        double sum = 0.0;
        for (int i = 0; i < ncol; ++i)
            sum += abs_sum(source1, source2, i);
        return sum;
    }

    double abs_sum(matrix<double> &source1, matrix<double> &source2, const int qi)
    {
        const int nrow = static_cast<int>(source1.n_rows);
        double sum = 0.0;
        for (int i = 0; i < nrow; ++i)
            sum += std::abs(source1(i, qi) - source2(i, qi));
        return sum;
    }

    double softnorm(const double s, const double lalph,
                    const double norm)
    {
        const double sabs = std::abs(s);
        if (lalph < sabs)
        {
            if (s > 0)
                return (s - lalph) / norm;
            return (s + lalph) / norm;
        }
        return 0.0;
    }

    double abs_dprod(const cvector<double> &lhs, const cvector<double> &rhs)
    {
        return std::abs(arma::dot(lhs, rhs));
    }

}