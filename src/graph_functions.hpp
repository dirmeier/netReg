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

#ifndef NETREG_GRAPH_FUNCTIONS_HPP
#define NETREG_GRAPH_FUNCTIONS_HPP

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

namespace netreg
{
    /**
     * Compute the degrees of the nodes of a graph.
     *
     * @param x a matrix
     *
     * @return the degrees of every node
     */
    std::vector<double> degree_distribution(const arma::Mat<double>& x);

    /**
     * Calculate the normalized Laplacian of a matrix.
     *
     * @param x a matrix
     *
     * @return the normalized laplacian
     */
    arma::Mat<double> laplacian(const arma::Mat<double>& x);
}
#endif  // NETREG_GRAPH_FUNCTIONS_HPP
