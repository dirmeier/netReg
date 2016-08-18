/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_MATH_HPP
#define NETREG_MATH_HPP

#include <vector>
#include <cmath>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "types.hpp"

namespace netreg
{
    /**
     * Calculate the dot product between two columns of two matrices and takes the absolute.
     *
     * @param source1 a random (n x p)-dimensional matrix
     * @param source2 a random (n x q)-dimensional matrix
     * @param pi the column of matrix<double> source1
     * @param qi the column of matrix<double> source2
     */
    double abs_dprod(const cvector<double> &lhs,
                     const cvector<double> &rhs);

    /**
     * Calculate a soft-thresholded normalized coefficient.
     *
     * @param s the estimate of the coefficient that should be estimated
     * @param lalph the Elastic-net regularization parameter
     * @param norm normalization constant
     * @return returns soft-thresholded normalized version of current coefficient
     */
    double softnorm(const double s,
                    const double lalph,
                    const double norm);

    /**
     * Returns the maximal element of a ptr
     *
     * @param ptr the pointer for which the maximum should be found
     * @param len length of the pointer
     *
     * @return return the maximal element
     */
    template <typename T> T max_element(T *const ptr, int len);


    /**
     * Calculates the sigmoid function value of a double.
     *
     * @param d  the value for which the sigmoidal is calculated
     *
     * @return returns the sigmoid function value
     */
    double sigmoid(double d)
    {
        return 1 / (1 + exp(d));
    }

}
#endif //NETREG_MATH_HPP
