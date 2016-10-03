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
#ifndef NETREG_GRAPH_PENALIZED_LINEAR_MODEL_CV_DATA_HPP
#define NETREG_GRAPH_PENALIZED_LINEAR_MODEL_CV_DATA_HPP

#include "graph_penalized_linear_model_data.hpp"

#include <vector>
#include <string>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "cv_set.hpp"
#include "types.hpp"

namespace netreg
{
    /**
     * Data-structure for all required data for graph-penalized regression.
     */
    class graph_penalized_linear_model_cv_data
        : public graph_penalized_linear_model_data
    {

    public:

        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param gx a pxp-dimensional prior graph for the covariables of x
         * @param gy a qxq-dimensional prior graph for the responses of y
         * @param n the number of samples (nrows X/Y)
         * @param p the number of covariables (ncols X)
         * @param q the number of responses (ncol Y)
         * @param lambda a vector of length q of penalisation values for q univariate models
         * @param alpha a vector of length q of weightings for lasso/ridge
         * @param psi_gx a vector of length q of how much influence GX should have on the penalization
         * @param psi_gy a vector of length q of how much influence GY should have on the penalization
         * @param niter max number of iterations in case estimation of the coefficients does not converge
         * @param thresh convergence threshold
         * @param nfolds the number of folds
         * @param fold_ids fold id mappings
         */
        graph_penalized_linear_model_cv_data
            (double *const x, double *const y,
             double *const gx, double *const gy,
             const int n, const int p, const int q,
             const double lambda, const double alpha,
             const double psi_gx, const double psi_gy,
             const int niter, const double thresh,
             const int nfolds, const enum family fam)
            : graph_penalized_linear_model_data
                  (x, y, gx, gy, n, p, q, lambda, alpha, psi_gx, psi_gy,
                   niter, thresh, fam),
              fold_ids_(design().n_rows), cvset_(n, nfolds)
        {
            set_fold_ids();
        }

        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param gx a pxp-dimensional prior graph for the covariables of x
         * @param gy a qxq-dimensional prior graph for the responses of y
         * @param n the number of samples (nrows X/Y)
         * @param p the number of covariables (ncols X)
         * @param q the number of responses (ncol Y)
         * @param lambda a vector of length q of penalisation values for q univariate models
         * @param alpha a vector of length q of weightings for lasso/ridge
         * @param psi_gx a vector of length q of how much influence GX should have on the penalization
         * @param psi_gy a vector of length q of how much influence GY should have on the penalization
         * @param niter max number of iterations in case estimation of the coefficients does not converge
         * @param thresh convergence threshold
         * @param fold_ids fold id mappings
         */
        graph_penalized_linear_model_cv_data
            (double *const x, double *const y,
             double *const gx, double *const gy,
             const int n, const int p, const int q,
             double const lambda, const double alpha,
             const double psi_gx, const double psi_gy,
             const int niter, const double thresh,
             int *const fold_ids, const enum family fam)
            : graph_penalized_linear_model_data
                  (x, y, gx, gy, n, p, q, lambda, alpha, psi_gx, psi_gy,
                   niter, thresh, fam),
              fold_ids_(design().n_rows), cvset_(n, fold_ids)
        {
            set_fold_ids();
        }

        /**
         * Getter for the vector of fold id mappings.
         *
         * @return returns the fold ids
         */
        std::vector<int> &fold_ids()
        {
            return fold_ids_;
        }

        /**
         * Getter
         */
        cv_set &cvset()
        {
            return cvset_;
        }

    private:
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

#endif //NETREG_GRAPH_PENALIZED_LINEAR_MODEL_CV_DATA_HPP
