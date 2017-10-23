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


#ifndef NETREG_GRAPHPENALIZEDLINEARMODELDATA_HPP
#define NETREG_GRAPHPENALIZEDLINEARMODELDATA_HPP

#include "penalized_linear_model_data.hpp"

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
    class graph_penalized_linear_model_data
        : public penalized_linear_model_data
    {

    public:
        graph_penalized_linear_model_data()
        {}

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
         */
        graph_penalized_linear_model_data
            (double *const x, double *const y,
             double *const gx, double *const gy,
             const int n, const int p, const int q,
             const double lambda, const double alpha,
             const double psi_gx, const double psi_gy,
             const int niter, const double thresh,
             const enum family fam)
            : penalized_linear_model_data(x, y, n, p, q, lambda,
                                          alpha, niter, thresh, fam),
              psi_gx(psi_gx), psi_gy(psi_gy),
              GX(gx, p, p, false, true), GY(gy, q, q, false, true),
              LX(laplacian(gx, p, p, psi_gx)),
              LY(laplacian(gy, q, q, psi_gy)),
              lx_rows_(LX.n_rows)
        {
            // copies LX rows as single vectors so that access can be faster
            #pragma omp parallel for
            for (std::vector<arma::Row<double> >::size_type i = 0;
                 i < LX.n_rows; ++i)
            {
                lx_rows_[i] = LX.row(i);
            }
        }

        /**
         * Getter for penalization term for laplacian of X
         *
         * @return the penalization term for X
         */
        const double psigx()
        {
            return psi_gx;
        }

        /**
         * Getter for penalization term for laplacian of Y
         *
         * @return the penalization term for Y
         */
        const double psigy()
        {
            return psi_gy;
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

        std::vector<arma::Row<double> >& lx_rows()
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
        double psi_gx;  // Penalization vector for GX
        double psi_gy;  // Penalization vector for GY
        arma::Mat<double> GX;    // prior graph for the design matrix
        arma::Mat<double> GY;    // prior graph for response matrix
        arma::Mat<double> LX;    // Normalized Laplacian of GX
        arma::Mat<double> LY;    // Normalized Laplacian of GY
        std::vector<arma::Row<double>> lx_rows_;
    };
}
#endif //NETREG_GRAPHPENALIZEDLINEARMODELDATA_HPP
