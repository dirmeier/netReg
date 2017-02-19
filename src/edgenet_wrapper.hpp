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

#include "graph_penalized_linear_model_data.hpp"
#include "graph_penalized_linear_model_cv_data.hpp"
#include "cv_fold.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>


namespace netreg
{
    class edgenet_wrapper
    {
    public:

        /**
         * Calulates the coefficients of a graph-regularized regression model.
         *
         * @param data an object that holds all required data for the model
         */
        SEXP run(graph_penalized_linear_model_data &data) const;

        /**
         * Calulates the optimal set of shrinkage parameters of a
         * graph-regularized regression model.
         *
         * @param data an object that holds all required data for the model
         * @param lambda the shrinkage parameter you want to use for the LASSO
         * @param alpha the parameter for the elastic net
         * @param psigx penalization of laplacian for X
         * @param psigy penalization of laplacian for Y
         *
         * @return returns the estimated coefficients
         */
       arma::Mat<double> run_cv
            (graph_penalized_linear_model_cv_data &data,
             const double lambda, const double alpha,
             const double psigx,  const double psigy,
             cv_fold &fold) const;
    };
}
#endif //NETREG_EDGENET_HPP
