/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "error_functions.hpp"

namespace netreg
{
    double sse(matrix<double> &B,
               matrix<double> &X,
               matrix<double> &Y,
               index_vector &testSet)
    {
        // sum of all elements of residuals
        double sum = arma::accu(Y.rows(testSet) - X.rows(testSet) * B);
        double sse = sum * sum;
        return sse;
    }
}