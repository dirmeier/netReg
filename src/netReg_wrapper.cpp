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
#include <memory>

#include "graph_penalized_linear_model_data.hpp"
#include "graph_penalized_linear_model_cv_data.hpp"
#include "edgenet_wrapper.hpp"
#include "edgenet_model_selection_wrapper.hpp"
#include "family.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
        Rprintf("Wrong family given\n");
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

extern "C"
{
/**
* Implementation of cross-validation for Edgenet.
*
* Finds and returns the optimal shrinkage values given a specific data-set.
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
* @param nfolds the number of cross-validation sets created (as in k-fold cv)
* @param foldids integer vector of assignments of observations
*        to folds (i.e. vector of ns elements,  \in {1, ..., nfolds}
* @param len_foldids length of the vector above
* @param fs family of distribution the response
*/
SEXP cv_edgenet_cpp
    (SEXP X, SEXP Y, SEXP GX, SEXP GY,
     SEXP psigx, SEXP psigy,
     SEXP niter, SEXP thresh,
     SEXP nfolds, SEXP foldids, SEXP lenfs,
     SEXP fs)
{
    BEGIN_RCPP;
    std::string fam = Rcpp::as<std::string>(fs);
    netreg::family f = fam == "gaussian" ? netreg::family::GAUSSIAN
                                         : netreg::family::NONE;
    if (f == netreg::family::NONE)
    {
        Rprintf("Wrong family given\n");
        return R_NilValue;
    }
    const int *xdim = INTEGER(Rf_getAttrib(X, R_DimSymbol));
    const int *ydim = INTEGER(Rf_getAttrib(Y, R_DimSymbol));
    const int lenfoldid = Rcpp::as<int>(lenfs);
    netreg::edgenet_model_selection e;
    if (lenfoldid == xdim[0])
    {
        netreg::graph_penalized_linear_model_cv_data data(
            REAL(X), REAL(Y), REAL(GX), REAL(GY), xdim[0], xdim[1], ydim[1],
            -1, 1.0, Rcpp::as<double>(psigx), Rcpp::as<double>(psigy),
            Rcpp::as<int>(niter), Rcpp::as<double>(thresh),
            INTEGER(foldids), f);
        return e.regularization_path(data);
    }
    netreg::graph_penalized_linear_model_cv_data data(
        REAL(X), REAL(Y), REAL(GX), REAL(GY), xdim[0], xdim[1], ydim[1],
        -1, 1.0, Rcpp::as<double>(psigx), Rcpp::as<double>(psigy),
        Rcpp::as<int>(niter), Rcpp::as<double>(thresh),
        Rcpp::as<int>(nfolds), f);
    return e.regularization_path(data);
    END_RCPP;
}
};


#include <R_ext/Rdynload.h>

static R_CallMethodDef callMethods[] = {
    {"edgenet_cpp", (DL_FUNC) &edgenet_cpp, 10},
    {"cv_edgenet_cpp", (DL_FUNC) &cv_edgenet_cpp, 12},
    {NULL, NULL, 0}
};


extern "C" void R_init_netReg(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
}
