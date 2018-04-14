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

#ifndef NETREG_EDGENET_HPP
#define NETREG_EDGENET_HPP

#include "graph_model_data.hpp"
#include "graph_model_cv_data.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace netreg
{
    /**
     * Calulates the coefficients of a graph-regularized regression model.
     *
     * @param data an object that holds all required data for the model
     * @param lamdba penalization value for LASSO
     * @param psigx weighting value of GX
     * @param psigy weighting value of GY
     * @param niter max number of iterations if parameter estimation
     *  does not converge in time
     * @param thresh convergence threshold
     */
    SEXP fit(
      graph_model_data& data,
      double lambda, double psigx, double psigy,
      int niter, double thresh);

    /**
     * Calulates the optimal set of shrinkage parameters of a
     * graph-regularized regression model.
     *
     * @param data an object that holds all required data for the model
     * @param lamdba penalization value for LASSO
     * @param psigx weighting value of GX
     * @param psigy weighting value of GY
     * @param do_lambda do estimation of lambda
     * @param do_psigx do estimation of psigx
     * @param do_psigy do estimation of psigy
     * @param niter max number of iterations if parameter estimation
     *        does not converge in time
     * @param thresh convergence threshold
     * @param optim_niter maximla number of iterations of BOBYQA
     * @param optim_epsilon threshold for convergence for BOBYQA
     *
     * @return returns the estimated parameters and foldids
     */
    SEXP regularization_path(
      graph_model_cv_data& data,
      double lambda, double psigx, double psigy,
      bool do_lambda, bool do_psigx, do_psigy,
      int niter, double thresh,
      int optim_niter, double optim_epsilon);
}
#endif  // NETREG_EDGENET_HPP
