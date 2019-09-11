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
     * Calculate the dot product between two matrices and takes
     * the absolute.
     *
     * @param lhs a random (n x p)-dimensional matrix
     * @param rhs a random (p x m)-dimensional matrix
     */
    inline double abs_dprod(const arma::Col<double>& lhs,
                            const arma::Col<double>& rhs)
    {
        return std::abs(arma::dot(lhs, rhs));
    }

    /**
     * Calculate a soft-thresholded normalized coefficient.
     *
     * @param s the value that will be soft-thresholded
     * @param lalph soft-thresholding parameter
     * @param norm normalization constant
     *
     * @return returns soft-thresholded value
     */
    inline double softnorm(const double s,
                           const double lalph,
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

    /**
     * Returns the maximal element of a ptr
     *
     * @param ptr the pointer for which the maximum should be found
     * @param len length of the pointer
     *
     * @return return the maximal element
     */
    template<typename T>
    T max_element(T* const ptr, int len)
    {
        T maximum = ptr[len - 1];
        for (unsigned int i = 0; i < len - 1; ++i)
        {
            if (ptr[i] > maximum)
                maximum = ptr[i];
        }

        return maximum;
    }

    /**
     * Calculates the sigmoid function.
     */
    inline double sigmoid(double d)
    {
        return 1 / (1 + exp(d));
    }

    /**
    * Computes the L1 norm of two vectors
    */
    inline double l1(const arma::vec& v1, const arma::vec& v2)
    {
        return arma::accu(arma::abs(v1 - v2));
    }

}
#endif  // NETREG_MATH_HPP
