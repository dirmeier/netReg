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
    class graph_model_cv_data
    {
    public:

        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param gx a pxp-dimensional prior graph for the covariables of x
         * @param gy a qxq-dimensional prior graph for the covariables of y
         * @param nfolds the number of folds
         * @param fam the fmily of the reponse variable
         *
         */
        graph_model_cv_data(
          arma::Mat<double>& x, arma::Mat<double>& y,
          arma::Mat<double>& gx, arma::Mat<double>& gy,
          int nfolds, const enum family fam):
            DATA_(x, y, gx, gy, fam),
            FOLD_IDS_(DATA_.sample_count()),
            CVSET_(DATA_.sample_count(), nfolds, DATA_.design(), DATA_.response())
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
         * @param fam the fmily of the reponse variable
         *
         */
        graph_model_cv_data(
          arma::Mat<double>& x, arma::Mat<double>& y,
          arma::Mat<double>& gx, arma::Mat<double>& gy,
          int nfolds, const enum family fam):
          DATA_(x, y, gx, gy, fam),
          FOLD_IDS_(DATA_.sample_count()),
          CVSET_(DATA_.sample_count(), nfolds, DATA_.design(), DATA_.response())
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
          enum family fam):
          DATA_(x, y, gx, gy, fam),
          FOLD_IDS_(DATA_.sample_count()),
          CVSET_(DATA_.sample_count(), fold_ids, DATA_.design(), DATA_.response())
        {
            throw not_implemented_exception();
            set_fold_ids();
        }

        /**
         * Getter for the vector of fold id mappings.
         */
        std::vector<int>& fold_ids()
        {
            return FOLD_IDS_;
        }

        /**
         * Getter for the graph model data.
         */
        graph_model_data& data()
        {
            return DATA_;
        }

        /**
         * Getter for the set of cross validation folds.
         */
        cv_set& cvset()
        {
            return CVSET_;
        }

    protected:

        /**
         * Function to set the fold id mappings to sample indexes.
         */
        void set_fold_ids();

    private:
        // the modelling data object
         graph_model_data DATA_;
        // mapping from fold id to index in samples
         std::vector<int> FOLD_IDS_;
        // the cross validation folds
         cv_set CVSET_;
    };
}

#endif
