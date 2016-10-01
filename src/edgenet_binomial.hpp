/**
 * Author: Simon Dirmeier
 * Date: 17/08/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_BINOMIAL_EDGENET_HPP
#define NETREG_BINOMIAL_EDGENET_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "edgenet.hpp"

#include "types.hpp"
#include "graph_penalized_linear_model_data.hpp"
#include "cv_fold.hpp"

namespace netreg
{
    /**
     * Class for estimating the coeffiecients of a edge-regularized linear regression model.
     */
    class edgenet_binomial :public edgenet
    {
    public:
        /**
         * Calulates the coefficients of a graph-regularized regression model.
         *
         * @param data an object that holds all required data for the model
         */
        void run(graph_penalized_linear_model_data &data) const;

        /**
         * Calulates the coefficients of a graph-regularized regression model.
         *
         * @param data an object that holds all required data for the model
         * @param lambda the shrinkage parameter you want to use for the LASSO
         * @param alpha the parameter for the elastic net
         * @param psigx penalization of laplacian for X
         * @param psigy penalization of laplacian for Y
         */
        matrix<double> run_cv
            (graph_penalized_linear_model_data &data,
             const double lambda,
             const double alpha,
             const double psigx,
             const double psigy,
             cv_fold &fold) const;
    };
}

#endif //NETREG_BINOMIAL_EDGENET_HPP
