/**
 * netReg: graph-regularized linear regression models.
 * <p>
 * Copyright (C) 2015 - 2016 Simon Dirmeier
 * <p>
 * This file is part of netReg.
 * <p>
 * netReg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * netReg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with netReg. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author: Simon Dirmeier
 * @email: simon.dirmeier@gmx.de
 */
#ifndef NETREG_MATH_FUNCTIONS_HPP
#define NETREG_MATH_FUNCTIONS_HPP

#include <vector>
#include <cmath>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

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
    double abs_dprod(const arma::Col<double> &lhs,
                     const arma::Col<double> &rhs);

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
    double sigmoid(double d);

}
#endif //NETREG_MATH_HPP
