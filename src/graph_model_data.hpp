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

#ifndef NETREG_GRAPHMODELDATA_HPP
#define NETREG_GRAPHMODELDATA_HPP

#include "model_data.hpp"

#include <string>
#include <vector>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "family.hpp"
#include "graph_functions.hpp"

namespace netreg
{
    /**
     * Data-structure for all required data for graph-penalized regression.
     */
    class graph_model_data: public model_data
    {
    public:

        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param gx a pxp-dimensional prior graph for the covariables of x
         * @param gy a qxq-dimensional prior graph for the responses of y
         * @param lambda a vector of length q of penalisation values for q
         *  univariate models
         * @param alpha a vector of length q of weightings for lasso/ridge
         * @param psi_gx a vector of length q of how much influence GX should
         *  have on the penalization
         * @param psi_gy a vector of length q of how much influence GY should
         *  have on the penalization
         *  @param do_psigx if true uses regularization of GX and psigx. Otherwise no
         *  regularization will be used.
         * @param do_psigy if true uses regularization of GY and psigy. Otherwise no
         *  regularization will be used.
         * @param niter max number of iterations in case estimation of the
         *  coefficients does not converge
         * @param thresh convergence threshold
         */
        graph_model_data(
          arma::Mat<double>& x, arma::Mat<double>& y,
          arma::Mat<double>& gx, arma::Mat<double>& gy,
          const enum family fam):
          model_data(x, y, fam),
          GX(gx.memptr(), gx.n_rows, gx.n_cols, false, true),
          GY(gy.memptr(), gy.n_rows, gy.n_cols, false, true),
          LX(laplacian(GX)),
          LY(laplacian(GY)),
          lx_rows_(LX.n_rows)
        {
            // stores a copy of LX rows as single vectors,
            // such that access will be faster
            for (std::vector<arma::Row<double>>::size_type i = 0;
                 i < LX.n_rows;
                 ++i)
            {
                lx_rows_[i] = LX.row(i);
            }
        }

        /**
         * Getter for affinity matrix of X
         *
         * @return reference to affinity matrix
         */
        arma::Mat<double>& gx()
        {
            return GX;
        }

        /**
         * Getter for affinity matrix of X
         *
         * @return reference to affinity matrix
         */
        arma::Mat<double>& gy()
        {
            return GY;
        }

        /**
         * Getter for laplacian matrix of X
         *
         * @return reference to laplacian matrix
         */
        arma::Mat<double>& lx()
        {
            return LX;
        }

        std::vector<arma::Row<double>>& lx_rows()
        {
            return lx_rows_;
        }

        /**
         * Getter for laplacian matrix of X
         *
         * @return reference to laplacian matrix
         */
        arma::Mat<double>& ly()
        {
            return LY;
        }

    protected:
        arma::Mat<double> GX;  // prior graph for the design matrix
        arma::Mat<double> GY;  // prior graph for response matrix
        arma::Mat<double> LX;  // Normalized Laplacian of GX
        arma::Mat<double> LY;  // Normalized Laplacian of GY
        std::vector<arma::Row<double>> lx_rows_;
    };
}
#endif
