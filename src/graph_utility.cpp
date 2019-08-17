/**
 * netReg: graph-regularized linear regression models.
 * <p>
 * Copyright (C) 2015 - 2019 Simon Dirmeier
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


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
#include <cmath>
#include <vector>
#include <algorithm>


/' Get the weighted node degrees of a adjacency matrix.
//'
//' @noRd
//' @param x matrix for which node degrees are computed
//' @return returns the node degrees as vectors
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
std::vector<double> node_degrees_(const arma::Mat<double>& x)
{
    std::vector<double> degrees(x.n_rows);

    #pragma omp parallel for
    for (unsigned int i = 0; i < x.n_rows; ++i)
    {
        double rowSum = 0.0;
        for (unsigned j = 0; j < x.n_cols; ++j) rowSum += x(i, j);
        degrees[i] = rowSum;
    }

    return degrees;
}

/' Compute a normalized graph Laplacian
//'
//' @noRd
//' @param  x  matrix for which Laplacian is computed
//' @return  returns a matrix
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export]]
arma::Mat<double> laplacian_(const arma::Mat<double>& x)
{
    std::vector<double> degrees = node_degrees_(x);
    arma::Mat<double> lap(x.n_rows, x.n_cols);

    #pragma omp parallel for
    for (unsigned int i = 0; i < x.n_rows; ++i)
    {
        for (unsigned int j = 0; j < x.n_cols; ++j)
        {
            if (i == j && degrees[i] != 0)
                lap(i, j) = 1 - (x(i, j) / degrees[i]);
            else if (i != j && x(i, j) != 0)
                lap(i, j) = -x(i, j) / std::sqrt(degrees[i] * degrees[j]);
            else
                lap(i, j) = 0.0;
        }
    }

    return lap;
}