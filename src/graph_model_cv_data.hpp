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

#ifndef NETREG_GRAPH_MODEL_CV_DATA_HPP
#define NETREG_GRAPH_MODEL_CV_DATA_HPP

#include "graph_model_data.hpp"

#include <vector>
#include <string>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

#include "cv_set.hpp"
#include "not_implemented_exception.hpp"

namespace netreg
{
    /**
     * Data-structure for all required data for graph-penalized regression.
     */
    class graph_model_cv_data: public graph_model_data
    {
    public:

        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param gx a pxp-dimensional prior graph for the covariables of x
         * @param nfolds the number of folds
         * @param fam the fmily of the reponse variable
         *
         */
        graph_model_cv_data(
          arma::Mat<double>& x, arma::Mat<double>& y,
          arma::Mat<double>& gx, arma::Mat<double>& gy,
          int nfolds, const enum family fam):
            graph_model_data(x, y, gx, gy, fam),
            fold_ids_(N), cvset_(N, nfolds, X, Y)
        {
            set_fold_ids();
        }

        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param gx a pxp-dimensional prior graph for the covariables of x
         * @param nfolds the number of folds
         * @param fold_ids fold id mappings
         * @param fam the fmily of the reponse variable
         *
         */
        graph_model_cv_data(
          arma::Mat<double>& x, arma::Mat<double>& y,
          arma::Mat<double>& gx, arma::Mat<double>& gy,
          int nfolds, int* const fold_ids,
          const enum family fam):
          graph_model_data(x, y, gx, gy, fam),
          fold_ids_(N), cvset_(N, fold_ids, X, Y)
        {
            throw not_implemented_exception();
            set_fold_ids();
        }

        /**
         * Getter for the vector of fold id mappings.
         */
        std::vector<int>& fold_ids()
        {
            return fold_ids_;
        }

        /**
         * Getter for the set of cross validation folds.
         */
        cv_set& cvset()
        {
            return cvset_;
        }

    protected:

        /**
         * Function to set the fold id mappings to sample indexes.
         */
        void set_fold_ids();
        // mapping from fold id to index in samples
        std::vector<int> fold_ids_;
        // the cross validation folds
        cv_set cvset_;
    };
}

#endif
