/**
 * Author: Simon Dirmeier
 * Date: 9/4/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_EDGENET_HPP
#define NETREG_EDGENET_HPP

#include "graph_penalized_linear_model_data.hpp"
#include "graph_penalized_linear_model_cv_data.hpp"
#include "cv_fold.hpp"

namespace netreg
{
    class edgenet
    {
    public:

        /**
         * Calulates the coefficients of a graph-regularized regression model.
         *
         * @param data an object that holds all required data for the model
         */
        virtual void run(graph_penalized_linear_model_data &data) const;

        /**
         * Calulates the optimal set of shrinkage parameters of a
         * graph-regularized regression model.
         *
         * @param data an object that holds all required data for the model
         * @param lambda the shrinkage parameter you want to use for the LASSO
         * @param alpha the parameter for the elastic net
         * @param psigx penalization of laplacian for X
         * @param psigy penalization of laplacian for Y
         *
         * @return returns the estimated coefficients
         */
       virtual matrix<double> run_cv
            (graph_penalized_linear_model_cv_data &data,
             const double lambda, const double alpha,
             const double psigx,  const double psigy,
             cv_fold &fold) const;
    };
}
#endif //NETREG_EDGENET_HPP
