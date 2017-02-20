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

#include "stat_functions.hpp"

namespace netreg
{
    arma::Col<double> intercept(arma::Mat<double> &X,
                                arma::Mat<double> &Y,
                              arma::Mat<double> &B)
    {
        arma::Mat<double> terr = (Y - (X * B)).t();
        arma::Col<double> rep(Y.n_rows, arma::fill::ones);
        arma::Col<double> intr = (terr * rep) / Y.n_rows;
        return intr;
    }

    double pls(arma::Mat<double> &TXX,
               arma::Mat<double> &TXY,
               arma::Mat<double> &cfs,
               const int pi,
               const int qi,
               const int P,
               std::vector< arma::rowvec >& txx_rows)
    {
            return TXY(pi, qi) + (TXX(pi, pi) * cfs(pi, qi))
                       - arma::accu(txx_rows[pi] * cfs.col(qi));
            
    }
}
