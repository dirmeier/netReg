/**
 * Author: Simon Dirmeier
 * Date: 17/08/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include "edgenet.hpp"
#include "math_functions.hpp"
#include "stat_functions.hpp"
#include <iostream>

namespace netreg
{
    void binomial_edgenet::run(graph_penalized_linear_model_data &data) const
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
        std::cout << "Running binomial edgenet\n";

        // calculate intercepts of the linear model
        std::cout << "Calculating intercept!\n";
        intr = intercept(data.design(), data.response(), coef);
        std::cout << "Done!\n";
    }

    matrix<double> binomial_edgenet::run_cv(graph_penalized_linear_model_data &data,
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

}