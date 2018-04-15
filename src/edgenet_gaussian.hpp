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

#ifndef NETREG_EDGENET_GAUSSIAN_HPP
#define NETREG_EDGENET_GAUSSIAN_HPP

#include <vector>
#include "params.hpp"
#include "graph_model_data.hpp"
#include "graph_model_cv_data.hpp"
#include "cv_fold.hpp"
#include "math_functions.hpp"
#include "stat_functions.hpp"

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

namespace netreg
{
    /**
     * Class for estimating the coefficients of a edge-regularized linear
     * regression model.
     */
    class edgenet_gaussian
    {
    public:

        edgenet_gaussian(const graph_model_data& data, const params& pars):
          DATA_(data),
          LX_(data.lx_rows()),
          LY_(data.ly()),
          P_(data.covariable_count()),
          Q_(data.response_count()),
          lambda_(pars.lambda()),
          psigx_(pars.psigx()),
          psigy_(pars.psigy()),
          NITER_(pars.niter()),
          THRESH_(pars.thresh())

        /**
         * Calulates the coefficients of a graph-regularized regression model.
         *
         * @param data an object that holds all required data for the model
         */
        arma::Mat<double> run() const;

        /**
         * Calculates the optimal set of shrinkage parameters of a
         * graph-regularized regression model.
         *
         * @param lambda the shrinkage parameter you want to use for the LASSO
         * @param psigx penalization of laplacian for X
         * @param psigy penalization of laplacian for Y
         * @param fold the fold you want to use of the data
         *
         * @return returns the estimated coefficients
         */
        arma::Mat<double> run_cv(double lambda,
                                 double psigx,
                                 double psigy,
                                 cv_fold& fold) const;


        void set_lambda(double lambda)
        {
            lambda_ = lambda;
        }

        void set_psigx()
        {
            psigx_ = psigx;
        }

        void set_psigy()
        {
            psigy = psigy;
        }

    protected:
        arma::Mat<double> mccd_(std::vector<arma::rowvec>& txx_rows,
                                arma::Mat<double>& txy) const;

        void uccd_(int qi,
                   std::vector<arma::rowvec>& txx_rows,
                   arma::Mat<double>& txy,
                   arma::Mat<double>& B,
                   arma::Mat<double>& B_old,
                   std::vector<arma::rowvec>& B_rows) const;

        inline void set_params(double& s,
                               double& norm,
                               const int pi,
                               const int qi,
                               arma::rowvec& txx_row,
                               arma::Mat<double>& txy,
                               arma::Mat<double>& B,
                               arma::rowvec& B_row) const
        {
            s = partial_least_squares(txx_row, txy, B, pi, qi);
            norm = txx_row(pi);
            graph_penalize(s, norm, pi, qi, B, B_row);
        }

        inline void graph_penalize(double& s,
                                   double& norm,
                                   const int pi,
                                   const int qi,
                                   arma::Mat<double>& B,
                                   arma::rowvec& B_row) const
        {
            if (psigx_ > 0.001 && LX_.size() == P_)
                lx_penalize(s, norm, pi, qi, psigx, B);
            if (psigy_ > 0.001 && LY_.n_rows == Q_ && Q_ > 1)
                ly_penalize(s, norm, pi, qi, psigy, B_row);
        }

        /**
         * Adds the penalty from the Laplacian of X to the softthresholding
         * parameter s and the normalization constant norm.
         *
         * @param s the softthresholding parameter to be set
         * @param norm the normalization constant to be set
         * @param cfs the current estimate of the coefficients
         * @param LX the normalized laplacian of X
         * @param P the number of covariables
         * @param pi the current index of the column of X
         * @param qi the current index of the column of Y
         * @param psigx the penalty for the Laplacian of X
         */
        inline void lx_penalize(double& s,
                                double& norm,
                                const int pi,
                                const int qi,
                                arma::Mat<double>& B) const
        {
            double xPenalty =
              -LX_[pi](pi) * B(pi, qi) + arma::accu(LX_[pi] * cfs.col(qi));
            s = s - 2 * psigx_ * xPenalty;
            norm += 2 * psigx_ * LX_[pi](pi);
        }

        /**
         * Adds the penalty from the Laplacian of X to the softthresholding
         * parameter s and the normalization constant norm.
         *
         * @param s the softthresholding parameter to be set
         * @param norm the normalization constant to be set
         * @param psigy the penalty for the Laplacien of Y
         * @param LY the normalized laplacian of Y
         * @param B the current estimate of the coefficients
         * @param Q the number of responses
         * @param pi the current index of the column of X
         * @param qi the current index of the column of Y
         */
        inline void ly_penalize(double& s,
                                double& norm,
                                const int pi,
                                const int qi,
                                arma::rowvec& B_row) const
        {
            double yPenalty = -B_row(qi) * LY_(qi, qi)
                              + arma::accu(B_row * LY_.col(qi));
            s = s - 2 * psigy * yPenalty;
            norm += 2 * psigy * LY_(qi, qi);
        }

    private:
        const graph_model_data& DATA_;
        const std::vector<arma::Row<double>>& LX_;
        const arma::Mat<double>& LY_;
        const int P_;
        const int Q_;
        double lambda_;
        double psigx_;
        double psigy_;
        const int NITER_;
        const double THRESH_;

    };
}

#endif  // NETREG_EDGENET_H
