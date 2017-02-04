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
#ifndef NETREG_EDGENET_MODEL_SELECTION_HPP
#define NETREG_EDGENET_MODEL_SELECTION_HPP

#include "graph_penalized_linear_model_cv_data.hpp"
#include "cv_set.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace netreg
{
    /**
     * Class that selects the best parameters of a model-selection using cross-validation
     */
    class edgenet_model_selection
    {
    public:

        /**
         * Calculate the optimal shrinkage parameters for an
         * edge-penalized regression model, i.e. calculate different
         * lambdas and psis.
         *
         * @param data the model data for which you want to estimate the
         *        optimal regularization parameters.
         *
         * @returns returns a SEXP with the optimal shrinkage parameters
         */
        SEXP regularization_path
            (graph_penalized_linear_model_cv_data &data);
    };
}
#endif //NETREG_EDGENETMODELSELECTION_HPP
