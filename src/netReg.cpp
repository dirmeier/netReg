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

#include <string>
#include "graph_penalized_linear_model_data.hpp"
#include "edgenet.hpp"
#include "family.hpp"

#include "Rcpp.h"

extern "C"
{
/**
 * Implementation of Edgenet, a edge-based regularized regression model.
 *
 * @param X a (n x p)-dimensional design matrix
 * @param Y a (n x q)-dimensional response matrix
 * @param GX a (p x p)-prior graph for XS
 * @param GY a (q x q)-prior graph for YS
 * @param lamdba penalization value for LASSO
 * @param psigx weighting value of GX
 * @param psigy weighting value of GY
 * @param niter max number of iterations if parameter estimation
 *        does not converge in time
 * @param thresh convergence threshold
 * @param fs family of distribution the response
 */
SEXP edgenet_cpp
    (SEXP X, SEXP Y, SEXP GX, SEXP GY,
     SEXP lambda, SEXP psigx, SEXP psigy,
     SEXP niter, SEXP thresh, SEXP fs)
{
    BEGIN_RCPP;
    std::string fam = Rcpp::as<std::string>(fs);
    netreg::family f = //fam == "binomial" ? netreg::family::BINOMIAL :
                       fam == "gaussian" ? netreg::family::GAUSSIAN
                                         : netreg::family::NONE;
    if (f == netreg::family::NONE)
    {
        Rcpp::Rcerr << "Wrong family given!" << "\n";
        return R_NilValue;
    }
    const int *xdim = INTEGER(Rf_getAttrib(X, R_DimSymbol));
    const int *ydim = INTEGER(Rf_getAttrib(Y, R_DimSymbol));
    netreg::graph_penalized_linear_model_data data
        (REAL(X), REAL(Y), REAL(GX), REAL(GY),
         xdim[0], xdim[1], ydim[1],
         Rcpp::as<double>(lambda), 1.0,
         Rcpp::as<double>(psigx), Rcpp::as<double>(psigy),
         Rcpp::as<int>(niter), Rcpp::as<double>(thresh), f);
    // TODO change that back and include family in data
    netreg::edgenet edge;
//    return R_NilValue;
    return edge.run(data);
    END_RCPP;
}
};
