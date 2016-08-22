/**
 * Author: Simon Dirmeier
 * Date: 17/08/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_BINOMIAL_EDGENET_LIKELIHOOD_FUNCTION_HPP
#define NETREG_BINOMIAL_EDGENET_LIKELIHOOD_FUNCTION_HPP

#include <numeric>
#include <cmath>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#include <omp.h>

#include "types.hpp"
#include "graph_penalized_linear_model_data.hpp"
#include "../inst/dlib/matrix.h"

namespace netreg
{
    /**
     * Functor class representing the objective function of a edge-regularized regression model.
     */
    class binomial_edgenet_likelihood_function
    {
    public:
        /**
         * Creates an objective function object that can be used for minimization using dlib.
         *
         * @param data the complete dataset required for edge-regularized regression
         * @param cvset a cross-validation set
         */
        binomial_edgenet_likelihood_function
            (graph_penalized_linear_model_data &data) :
            data_(data),
            X_(data.design()),
            Y_(data.response()),
            nfolds_(static_cast<int>(cvset.fold_count())),
            edgenet_(),
            do_psigx_(data.psigx() == -1),
            do_psigy_(data.psigy() == -1),
            P_(X_.n_cols),
            Q_(Y_.n_cols)
        {}

        /**
         * Over-write operator () in order to get functor functionality (object behaves like a function)
         *
         * @param params the free parameters on the objective function
         */
        double operator()(const dlib::matrix<double> &b) const
        {

            matrix<double> B(P, Q);
            // TODO: efficient cast b to B
            for (unsigned int i = 0; i < P; ++i)
                for (unsigned int j = 0; j < Q; ++j)
                    B(i, j) = b(i, j);
            matrix<double> sigm = (X_ * B);
            sigm.transform([](double val) { return log(1 / (1 + exp(val))); });
            matrix<double> loglik = Y .* log(sigm)  + (1 - Y) .* log(1-sigm);
            double nll = 0.0;

            return nll;
        }

    private:
        // data required for a edge-regularized regression model
        graph_penalized_linear_model_data &data_;
        matrix<double> &X_;      // design matrix
        matrix<double> &Y_;      // response matrix
        const bool do_psigx_;
        const bool do_psigy_;
        const int P_;
        const int Q_;
    };
}

#endif //NETREG_BINOMIAL_EDGENET_LIKELIHOOD_FUNCTION_HPP
