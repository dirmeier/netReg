/**
 * Author: Simon Dirmeier
 * Date: 14/05/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_ERROR_FUNCTIONS_HPP
#define NETREG_ERROR_FUNCTIONS_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#include "../struc/types.hpp"

namespace netreg
{
    /**
     * Calculate the sum of squaerd errors.
     *
     * @param B the estimator for the coefficient matrix
     * @param X the design matrix
     * @param Y the response matrix
     * @param test_set a vector of indexed of the test set
     */
    double sse(matrix<double> &B,
               matrix<double> &X,
               matrix<double> &Y,
               std::vector<int> &test_set);
}
#endif //SRC_ERROR_FUNCTIONS_HPP
