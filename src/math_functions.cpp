/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "math_functions.hpp"

#include <cmath>

namespace netreg
{
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

    template<typename T>
    T max_element(T *const ptr, int len)
    {
        T maximum = ptr[len - 1];
        for (int i = 0; i < len - 1; ++i)
            if (ptr[i] > maximum) maximum = ptr[i];
        return maximum;
    }

    double sigmoid(double d)
    {
        return 1 / (1 + exp(d));
    }
}