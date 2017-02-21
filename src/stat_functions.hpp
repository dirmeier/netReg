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
#ifndef NETREG_STAT_FUNCTIONS_HPP
#define NETREG_STAT_FUNCTIONS_HPP

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

#include <vector>

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
    arma::Col<double> intercept(arma::Mat<double> &X,
                                arma::Mat<double> &Y,
                                arma::Mat<double> &B);

    /**
     * Calculates the partial residual of the current coefficient that is estimated.
     *
     * @param TXX the square of the design matrix
     * @param TXY the design times the response matrix
     * @param cfs the current estimate of the coefficients
     * @param pi the current index of the column of X
     * @param qi the current index of the column of Y
     */
    inline double partial_least_squares(arma::rowvec &txx_rows,
                                        arma::Mat<double> &TXY,
                                        arma::Mat<double> &cfs,
                                        const int pi, const int qi)
    {
        return TXY(pi, qi) + (txx_rows(pi) * cfs(pi, qi))
               - arma::accu(txx_rows * cfs.col(qi));
    }
}
#endif //NETREG_STAT_FUNCTIONS_HPP
