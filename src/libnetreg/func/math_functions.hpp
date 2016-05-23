#ifndef NETREG_MATH_HPP
#define NETREG_MATH_HPP

#include <vector>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "../struc/types.hpp"

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
     * Calculate the sum of the absolute differences between two matrices
     *
     * @param m1 a matrix object
     * @param m2 a matrix object
     * @return returns the sum of absolute differences
     */
    double abs_sum(matrix<double> &m1, matrix<double> &m2);

    /**
     * Calculate the sum of the absolute differences between two columns of two matrices
     *
     * @param m1 a matrix object
     * @param m2 a matrix object
     * @param col_idx the index of the columns
     * @return returns the sum of absolute differences
     */
    double abs_sum(matrix<double> &source1, matrix<double> &source2,
                   const int col_idx);

}
#endif //NETREG_MATH_HPP
