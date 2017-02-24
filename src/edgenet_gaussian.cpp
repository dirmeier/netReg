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

#include "edgenet_gaussian.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "cv_set.hpp"

namespace netreg
{
    arma::Mat<double> edgenet_gaussian::run(
        graph_penalized_linear_model_data &data) const
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();

        const double thresh = data.threshold();
        const int niter = data.max_iter();

        const double lambda = data.lambda();
        const double alpha = data.alpha();
        const double psigx = data.psigx();
        const double psigy = data.psigy();

        arma::Mat<double> &TXY = data.txy();
        arma::Mat<double> &LY = data.ly();

        std::vector <arma::rowvec> &txx_rows = data.txx_rows();
        std::vector <arma::Row<double>> &lx_rows = data.lx_rows();

        arma::Mat<double> coef(P, Q, arma::fill::ones);
        arma::Mat<double> old_coef(P, Q);
        std::vector <arma::rowvec> coef_rows(static_cast<unsigned int>(P));

        #pragma omp parallel for
        for (std::vector<arma::Row < double >>::size_type i = 0; i < coef.n_rows; ++i) {
            coef_rows[i] = coef.row(i);
        }

        // weighted penalization param of Elastic-net
        const double lalph = alpha * lambda;
        // normalization for soft-thresholding
        const double enorm = 1.0 + lambda * (1 - alpha);

        int iter = 0;
        do
        {
            // This is definitely not parallizable
            // do not try it in the future
            for (int qi = 0; qi < Q; ++qi)
            {
                uccd_(P, Q, qi,
                      thresh, niter,
                      lalph, enorm,
                      psigx, psigy,
                      TXY, LY,
                      coef, old_coef,
                      txx_rows,
                      lx_rows,
                      coef_rows);
#ifdef USE_RCPPARMADILLO
                if (iter % 100 == 0) Rcpp::checkUserInterrupt();
#endif
            }
        }
        while (arma::accu(arma::abs(coef - old_coef)) > thresh &&
               iter++ < niter);
        return coef;
    }

    arma::Mat<double> edgenet_gaussian::run_cv
        (graph_penalized_linear_model_cv_data &data,
         const double lambda,
         const double alpha,
         const double psigx,
         const double psigy,
         cv_fold &fold) const
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();

        const double thresh = data.threshold();
        const int niter = data.max_iter();

        arma::Mat<double> &X = data.design();
        arma::Mat<double> &Y = data.response();
        arma::Mat<double> &LY = data.ly();
        arma::uvec &trainIdxs = fold.train_set();

        arma::Mat<double> Xtrain = X.rows(trainIdxs);
        arma::Mat<double> Ytrain = Y.rows(trainIdxs);
        arma::Mat<double> TXtrain = Xtrain.t();
        arma::Mat<double> train_txx = TXtrain * Xtrain;
        arma::Mat<double> train_txy = TXtrain * Ytrain;

        std::vector <arma::rowvec> txx_rows(P);
        #pragma omp parallel for
        for (std::vector < arma::Row < double > > ::size_type i = 0; i < train_txx.n_rows; ++i) {
            txx_rows[i] = train_txx.row(i);
        }

        std::vector <arma::Row<double>> &lx_rows = data.lx_rows();

        arma::Mat<double> coef(P, Q, arma::fill::ones);
        arma::Mat<double> old_coef(P, Q);
        std::vector <arma::rowvec> coef_rows(P);
        #pragma omp parallel for
        for (std::vector < arma::Row < double > > ::size_type i = 0; i < coef.n_rows; ++i) {
            coef_rows[i] = coef.row(i);
        }

        // weighted penalization param of Elastic-net
        const double lalph = alpha * lambda;
        // normalization for soft-thresholding
        const double enorm = 1.0 + lambda * (1 - alpha);

        int iter = 0;
        do
        {
            // This is definitely not parallizable!
            for (int qi = 0; qi < Q; ++qi)
            {
                uccd_(P, Q, qi,
                      thresh, niter,
                      lalph, enorm,
                      psigx, psigy,
                      train_txy, LY,
                      coef, old_coef,
                      txx_rows,
                      lx_rows,
                      coef_rows);
#ifdef USE_RCPPARMADILLO
                if (iter % 100 == 0) Rcpp::checkUserInterrupt();
#endif
            }
        }
        while (arma::accu(arma::abs(coef - old_coef)) > thresh &&
               iter++ < niter);
        return coef;
    }

    void edgenet_gaussian::uccd_
        (const int P, const int Q, const int qi,
         const double thresh, const int niter,
         const double lalph, const double enorm,
         const double psigx, const double psigy,
         arma::Mat<double> &TXY,
         arma::Mat<double> &LY,
         arma::Mat<double> &coef,
         arma::Mat<double> &old_coef,
         std::vector <arma::rowvec> &txx_rows,
         std::vector <arma::rowvec> &lx_rows,
         std::vector <arma::rowvec> &coef_rows) const
    {
        // iteration counter
        int iter = 0;
        double s = 0.0;
        double norm = 0.0;
        // do while estimation of params does not converge
        do
        {
            // fix bnew_i and calculate least-squares
            // coefficient on partial residual
            for (int pi = 0; pi < P; ++pi)
            {
                // safe current estimate of coefficients
                old_coef(pi, qi) = coef(pi, qi);
                set_params
                    (s, norm, P, Q,
                     pi, qi, psigx, psigy,
                     TXY, LY, coef,
                     txx_rows[pi],
                     lx_rows[pi],
                     coef_rows[pi]);
                // soft-thresholded version of estimate
                const double d = softnorm(s, lalph, enorm * norm);
                coef(pi, qi) = d;
                coef_rows[pi](qi) = d;
#ifdef USE_RCPPARMADILLO
                if (iter % 100 == 0) Rcpp::checkUserInterrupt();
#endif
            }
        }
        while (
            arma::accu(arma::abs(coef.col(qi) - old_coef.col(qi))) > thresh &&
            iter++ < niter);
    }
}

