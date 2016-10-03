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

#ifndef NETREG_ERROR_FUNCTIONS_HPP
#define NETREG_ERROR_FUNCTIONS_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#include "types.hpp"

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
               index_vector &test_set);
}
#endif //SRC_ERROR_FUNCTIONS_HPP
