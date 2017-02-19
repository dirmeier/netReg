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

#include "cv_set.hpp"
#include "math_functions.hpp"
#include "stat_functions.hpp"

namespace netreg
{
    arma::Mat<double> edgenet_gaussian::run(graph_penalized_linear_model_data &data) const
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();
        arma::Mat<double> coef(P, Q, arma::fill::ones);
        arma::Mat<double> old_coef(P, Q);
        const double thresh = data.threshold();
        const int niter = data.max_iter();
        const double lambda = data.lambda();
        const double alpha = data.alpha();
        const double psigx = data.psigx();
        const double psigy = data.psigy();
        arma::Mat<double> &TXX = data.txx();
        arma::Mat<double> &TXY = data.txy();
        arma::Mat<double> &LX = data.lx();
        arma::Mat<double> &LY = data.ly();
        int iter = 0;
        do
        {
            // This is definitely not parallizable
            for (int qi = 0; qi < Q; ++qi)
            {
                uccd_(P, Q,
                      thresh, niter,
                      lambda, alpha,
                      psigx, psigy,
                      TXX, TXY,
                      LX, LY,
                      coef, old_coef,
                      qi);
#ifdef USE_RCPPARMADILLO
                if (iter % 10 == 0) Rcpp::checkUserInterrupt();
#endif
            }
        }
        while (arma::accu(arma::abs(coef - old_coef)) > thresh &&
               iter++ < niter);
        return coef;
    }

    arma::Mat<double> edgenet_gaussian::run_cv(
        graph_penalized_linear_model_cv_data &data,
        const double lambda, const double alpha,
        const double psigx, const double psigy,
        cv_fold &fold) const
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();
        arma::Mat<double> coef(P, Q, arma::fill::ones);
        arma::Mat<double> old_coef(P, Q);
        const double thresh = data.threshold();
        const int niter = data.max_iter();
        arma::Mat<double> &X = data.design();
        arma::Mat<double> &Y = data.response();
        arma::Mat<double> &LX = data.lx();
        arma::Mat<double> &LY = data.ly();
        arma::uvec &trainIdxs = fold.train_set();
        arma::Mat<double> Xtrain = X.rows(trainIdxs);
        arma::Mat<double> Ytrain = Y.rows(trainIdxs);
        arma::Mat<double> TXtrain = Xtrain.t();
        arma::Mat<double> train_txx = TXtrain * Xtrain;
        arma::Mat<double> train_txy = TXtrain * Ytrain;
        int iter = 0;
        do
        {
            // This is definitely not parallizable!
            for (int qi = 0; qi < Q; ++qi)
            {
                uccd_(P, Q,
                      thresh, niter,
                      lambda, alpha,
                      psigx, psigy,
                      train_txx, train_txy,
                      LX, LY,
                      coef, old_coef,
                      qi);
#ifdef USE_RCPPARMADILLO
                if (iter % 10 == 0) Rcpp::checkUserInterrupt();
#endif
            }
        }
        while (arma::accu(arma::abs(coef - old_coef)) > thresh &&
               iter++ < niter);
        return coef;
    }

    void edgenet_gaussian::uccd_
        (const int P, const int Q,
         const double thresh, const int niter,
         const double lambda, const double alpha,
         const double psigx, const double psigy,
         arma::Mat<double> &TXX, arma::Mat<double> &TXY,
         arma::Mat<double> &LX, arma::Mat<double> &LY,
         arma::Mat<double> &coef,
         arma::Mat<double> &old_coef,
         const int qi) const
    {
        // weighted penalization param of Elastic-net
        const double lalph = alpha * lambda;
        // normalization for soft-thresholding
        const double enorm = 1.0 + lambda * (1 - alpha);
        // iteration counter
        int iter = 0;
        // do while estimation of params does not converge
        do
        {
            // fix bnew_i and calculate least-squares
            // coefficient on partial residual
            for (int pi = 0; pi < P; ++pi)
            {
                // safe current estimate of coefficients
                old_coef(pi, qi) = coef(pi, qi);
                double s = 0.0;
                double norm = 0.0;
                set_params
                    (s, norm, TXX, TXY, coef,
                     LX, LY, P, Q, pi, qi, psigx, psigy, false);
//                // soft-thresholded version of estimate
                coef(pi, qi) = softnorm(s, lalph, enorm * norm);
#ifdef USE_RCPPARMADILLO
                if (iter % 10 == 0) Rcpp::checkUserInterrupt();
#endif
            }
        }
        while (
            arma::accu(arma::abs(coef.col(qi) - old_coef.col(qi))) > thresh &&
            iter++ < niter);
    }

    void edgenet_gaussian::set_params
        (double &s, double &norm,
         arma::Mat<double> &TXX, arma::Mat<double> &TXY,
         arma::Mat<double> &coef,
         arma::Mat<double> &LX,
         arma::Mat<double> &LY,
         const int P, const int Q,
         const int pi, const int qi,
         const double psigx,
         const double psigy,
         const bool lower) const
    {
        s = pls(TXX, TXY, coef, pi, qi, P, lower);
        norm = (TXX)(pi, pi);
        graph_penalize(s, norm, psigx, psigy,
                       LX, LY, coef,
                       P, Q, pi, qi);
    }

    void edgenet_gaussian::graph_penalize
        (double &s, double &norm,
         const double psigx, const double psigy,
         arma::Mat<double> &LX, arma::Mat<double> &LY, arma::Mat<double> &cfs,
         const int P, const int Q,
         const int pi, const int qi) const
    {
        if (psigx != 0)
            lx_penalize(s, norm, psigx, LX, cfs, P, pi, qi);
        if (psigy != 0)
            ly_penalize(s, norm, psigy, LY, cfs, Q, pi, qi);
    }

    void edgenet_gaussian::lx_penalize
        (double &s, double &norm, const double psigx,
         arma::Mat<double> &LX, arma::Mat<double> &cfs, const int P,
         const int pi, const int qi) const
    {
        if (psigx <= 0.001)
            return;
        double xPenalty = 0.0;
        if (pi < LX.n_rows && pi < LX.n_cols)
        {
            xPenalty = -LX(pi, pi) * cfs(pi, qi) + arma::accu(LX.row(pi) * cfs.col(qi));
        }
        s = s - 2 * psigx * xPenalty;
        norm += 2 * psigx * LX(pi, pi);
    }

    void edgenet_gaussian::ly_penalize
        (double &s, double &norm, const double psigy,
         arma::Mat<double> &LY, arma::Mat<double> &cfs, const int Q,
         const int pi, const int qi) const
    {
        if (psigy <= 0.001 || Q == 1)
            return;
        double yPenalty = 0.0;
        if (qi < LY.n_rows && qi < LY.n_cols)
        {
            yPenalty = -cfs(pi, qi) * LY(qi, qi) + arma::accu(cfs.row(pi) * LY.col(qi));
        }
        s = s - 2 * psigy * yPenalty;
        norm += 2 * psigy * LY(qi, qi);
    }

}
