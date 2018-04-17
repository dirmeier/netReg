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

#include "params.hpp"
#include "data_factory.hpp"
#include "edgenet_wrapper.hpp"
#include "graph_model_data.hpp"
#include "graph_model_cv_data.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using data = netreg::graph_model_data;
using cv_data = netreg::graph_model_cv_data;


extern "C" {

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
 *  regularization will be used.
 * @param niter max number of iterations if parameter estimation
 *  does not converge in time
 * @param thresh convergence threshold
 * @param fs family of distribution the response
 */
SEXP edgenet_cpp(SEXP X,
                 SEXP Y,
                 SEXP GX,
                 SEXP GY,
                 SEXP lambda, // 5
                 SEXP psigx,
                 SEXP psigy,
                 SEXP niter,
                 SEXP thresh,
                 SEXP fs)
{
    BEGIN_RCPP

    std::string f = Rcpp::as<std::string>(fs);
    data dat = netreg::data_factory::build_data(
      REAL(X), REAL(Y), REAL(GX), REAL(GY),
      INTEGER(Rf_getAttrib(X, R_DimSymbol)),
      INTEGER(Rf_getAttrib(Y, R_DimSymbol)),
      f
    );

    netreg::params pars = netreg::params()
      .lambda(Rcpp::as<double>(lambda))
      .psigx(Rcpp::as<double>(psigx))
      .psigy(Rcpp::as<double>(psigy))
      .thresh(Rcpp::as<double>(thresh))
      .niter(Rcpp::as<int>(niter));

    return fit(dat, pars);

    END_RCPP
    return R_NilValue;
}
};

extern "C" {

/**
 * Implementation of cross-validation for Edgenet.
 *
 * Finds and returns the optimal shrinkage values given a specific data-set.
 *
 * @param X a (n x p)-dimensional design matrix
 * @param Y a (n x q)-dimensional response matrix
 * @param GX a (p x p)-prior graph for XS
 * @param GY a (q x q)-prior graph for YS
 * @param lambda regularization parameter for LASSO
 * @param psigx weighting value of GX
 * @param psigy weighting value of GY
 * @param do_lambda do estimation of lambda
 * @param do_psigx do estimation of psigx
 * @param do_psigy do estimation of psigy
 * @param niter max number of iterations if parameter estimation
 *        does not converge in time
 * @param thresh convergence threshold
 * @param nfolds the number of cross-validation sets created (as in k-fold cv)
 * @param foldids integer vector of assignments of observations
 *        to folds (i.e. vector of ns elements,  \in {1, ..., nfolds}
 * @param lenfs length of the vector above
 * @param fs family of distribution the response
 * @param optim_niter maximla number of iterations of BOBYQA
 * @param optim_epsilon threshold for convergence for BOBYQA
 */
SEXP cv_edgenet_cpp(SEXP X,
                    SEXP Y,
                    SEXP GX,
                    SEXP GY,
                    SEXP lambda,  // 5
                    SEXP psigx,
                    SEXP psigy,
                    SEXP do_lambda,
                    SEXP do_psigx,
                    SEXP do_psigy, // 10
                    SEXP niter,
                    SEXP thresh,
                    SEXP nfolds,
                    SEXP foldids,
                    SEXP lenfs,  // 15
                    SEXP fs,
                    SEXP optim_niter,
                    SEXP optim_epsilon)
{
    BEGIN_RCPP

    std::string f = Rcpp::as<std::string>(fs);
    cv_data dat = netreg::data_factory::build_cv_data(
      REAL(X), REAL(Y), REAL(GX), REAL(GY),
      INTEGER(Rf_getAttrib(X, R_DimSymbol)),
      INTEGER(Rf_getAttrib(Y, R_DimSymbol)),
      f,
      Rcpp::as<int>(nfolds),
      Rcpp::as<int>(lenfs),
      INTEGER(foldids));

    netreg::params pars = netreg::params()
      .lambda(Rcpp::as<double>(lambda))
      .psigx(Rcpp::as<double>(psigx))
      .psigy(Rcpp::as<double>(psigy))
      .do_lambda(Rcpp::as<bool>(do_lambda))
      .do_psigx(Rcpp::as<bool>(do_psigx))
      .do_psigy(Rcpp::as<bool>(do_psigy))
      .thresh(Rcpp::as<double>(thresh))
      .niter(Rcpp::as<int>(niter))
      .optim_niter(Rcpp::as<int>(optim_niter))
      .optim_epsilon(Rcpp::as<double>(optim_epsilon));

    return regularization_path(dat, pars);

    END_RCPP

    return R_NilValue;
}
};

#include <R_ext/Rdynload.h>

static R_CallMethodDef callMethods[] = {
  {"edgenet_cpp",    (DL_FUNC) & edgenet_cpp,    10},
  {"cv_edgenet_cpp", (DL_FUNC) & cv_edgenet_cpp, 18},
  {NULL, NULL,                                   0}};

extern "C" void R_init_netReg(DllInfo* dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
