/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_STAT_FUNCTIONS_HPP
#define NETREG_STAT_FUNCTIONS_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "types.hpp"

namespace netreg
{
    /**
     * Calculate the intercept of a linear model
     *
     * @param X the design matrix
     * @param Y the response matrix
     * @param B the estimated coefficients
     * @return returns a column vector
     */
    cvector<double> intercept(matrix<double> &X,
                              matrix<double> &Y,
                              matrix<double> &B);

    /**
     * Calculate the partial residual with respect to pi given a row and a response.
     *
     * @param X the design matrix
     * @param Y the response matrix
     * @param B the estimated coefficient matrix
     * @param row the row for which the partial residual is calculated
     * @param pi the index of the coefficient
     * @param qi the index of the response
     */
    double partial_residual(matrix<double> &X,
                            matrix<double> &Y,
                            matrix<double> &B,
                            const int row,
                            const int pi,
                            const int qi);
    /**
     * Calculates the partial residual of the current coefficient that is estimated
     * IF ONLY THE LOWER TRIANGULAR MATRIX of TXX is given
     *
     * @param TXX the square of the design matrix
     * @param TXY the design times the response matrix
     * @param cfs the current estimate of the coefficients
     * @param P the number of covariables
     * @param pi the current index of the column of X
     * @param qi the current index of the column of Y
     */
    double l_pls(matrix<double> &TXX, matrix<double> &TXY,
                 matrix<double> &cfs,
                 const int cidx, const int qi, const int P);

    /**
     * Calculates the partial residual of the current coefficient that is estimated.
     *
     * @param TXX the square of the design matrix
     * @param TXY the design times the response matrix
     * @param cfs the current estimate of the coefficients
     * @param P the number of covariables
     * @param pi the current index of the column of X
     * @param qi the current index of the column of Y
     * @param lower is only lower of TXX is initialized
     */
    double pls(matrix<double> &TXX, matrix<double> &TXY, matrix<double> &cfs,
               const int pi, const int qi, const int P, const bool lower);
}
#endif //NETREG_STAT_FUNCTIONS_HPP
