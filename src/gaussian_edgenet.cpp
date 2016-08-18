/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "gaussian_edgenet.hpp"
#include "math_functions.hpp"
#include "stat_functions.hpp"
#include <iostream>

namespace netreg
{
    void gaussian_edgenet::run(graph_penalized_linear_model_data &data) const
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();
        matrix<double> &coef = data.coefficients();
        matrix<double> old_coef(P, Q);
        cvector<double> &intr = data.intercept();
        const double thresh = data.threshold();
        const int niter = data.max_iter();
        const double lambda = data.lambda();
        const double alpha = data.alpha();
        const double psigx=data.psigx();
        const double psigy=data.psigy();
        matrix<double> &TXX= data.txx();
        matrix<double> &TXY = data.txy();
        matrix<double> &LX = data.lx();
        matrix<double> &LY = data.ly();
        std::cout << "Running edgenet\n";
        int iter = 0;
        do
        {
            // This is definitely not parallizable
            for (int qi = 0; qi < Q; ++qi)
                uccd_(P, Q,
                      thresh, niter,
                      lambda, alpha,
                      psigx, psigy,
                      TXX, TXY,
                      LX, LY,
                      coef, old_coef,
                      qi);
        }
        while (arma::accu(arma::abs(coef - old_coef)) > thresh &&
               iter++ < niter);
        // calculate intercepts of the linear model
        std::cout << "Calculating intercept!\n";
        intr = intercept(data.design(), data.response(), coef);
        std::cout << "Done!\n";
    }

    matrix<double> gaussian_edgenet::run_cv(graph_penalized_linear_model_data &data,
                                  const double lambda, const double alpha,
                                  const double psigx, const double psigy,
                                  cv_fold &fold) const
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();
        matrix<double> coef(P, Q, arma::fill::ones);
        matrix<double> old_coef(P, Q);
        const double thresh = data.threshold();
        const int niter = data.max_iter();
        matrix<double> &X = data.design();
        matrix<double> &Y = data.response();
        matrix<double> &LX = data.lx();
        matrix<double> &LY = data.ly();
        index_vector &trainIdxs = fold.train_set();
        matrix<double> Xtrain = X.rows(trainIdxs);
        matrix<double> Ytrain = Y.rows(trainIdxs);
        matrix<double> TXtrain = Xtrain.t();
        matrix<double> train_txx = TXtrain * Xtrain;
        matrix<double> train_txy = TXtrain * Ytrain;
        int iter = 0;
        do
        {
            // This is definitely not parallizable!
            for (int qi = 0; qi < Q; ++qi)
            {
                uccd_(P, Q, thresh, niter, lambda, alpha, psigx, psigy,
                      train_txx, train_txy, LX, LY, coef, old_coef, qi);
            }
        }
        while (arma::accu(arma::abs(coef - old_coef)) > thresh && iter++ < niter);
        return coef;
    }

    void gaussian_edgenet::uccd_
        (const int P, const int Q,
         const double thresh, const int niter,
         const double lambda, const double alpha,
         const double psigx, const double psigy,
         matrix<double> &TXX, matrix<double> &TXY,
         matrix<double> &LX, matrix<double> &LY,
         matrix<double> &coef,
         matrix<double> &old_coef,
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
                // soft-thresholded version of estimate
                coef(pi, qi) = softnorm(s, lalph, enorm * norm);
            }
        }
        while (arma::accu(arma::abs(coef.col(qi) - old_coef.col(qi))) > thresh && iter++ < niter);
    }

    void gaussian_edgenet::set_params
        (double &s, double &norm,
         matrix<double> &TXX, matrix<double> &TXY,
         matrix<double> &coef,
         matrix<double> &LX,
         matrix<double> &LY,
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

    void gaussian_edgenet::graph_penalize
        (double &s, double &norm,
         const double psigx, const double psigy,
         matrix<double> &LX, matrix<double> &LY, matrix<double> &cfs,
         const int P, const int Q,
         const int pi, const int qi) const
    {
        if (psigx != 0)
            lx_penalize(s, norm, psigx, LX, cfs, P, pi, qi);
        if (psigy != 0)
            ly_penalize(s, norm, psigy, LY, cfs, Q, pi, qi);
    }

    void gaussian_edgenet::lx_penalize
        (double &s, double &norm, const double psigx,
         matrix<double> &LX, matrix<double> &cfs, const int P,
         const int pi, const int qi) const
    {
        if (psigx == 0)
            return;
        double xPenalty = -LX(pi, pi) * cfs(pi, qi);
        for (int j = 0; j < P; j++)
            xPenalty += LX(pi, j) * cfs(j, qi);
        s = s - 2 * psigx * xPenalty;
        norm += 2 * psigx * LX(pi, pi);
    }

    void gaussian_edgenet::ly_penalize
        (double &s, double &norm, const double psigy,
         matrix<double> &LY, matrix<double> &cfs, const int Q,
         const int pi, const int qi) const
    {
        if (psigy == 0 || Q == 1)
            return;
        // penalization for GY
        double yPenalty = -cfs(pi, qi) * LY(qi, qi);
        for (int j = 0; j < Q; ++j)
            yPenalty += cfs(pi, j) * LY(j, qi);
        s = s - 2 * psigy * yPenalty;
        norm += 2 * psigy * LY(qi, qi);
    }

}