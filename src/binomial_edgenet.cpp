/**
 * Author: Simon Dirmeier
 * Date: 17/08/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include "binomial_edgenet.hpp"
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

        matrix<double> coef(10, 10, arma::fill::ones);

        return coef;
    }

}