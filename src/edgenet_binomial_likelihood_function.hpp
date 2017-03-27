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
#ifndef NETREG_EDGENET_BINOMIAL_LIKELIHOOD_FUNCTION_HPP
#define NETREG_EDGENET_BINOMIAL_LIKELIHOOD_FUNCTION_HPP

#include <numeric>
#include <cmath>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "graph_penalized_linear_model_data.hpp"
#include "../inst/dlib/matrix.h"

namespace netreg
{
    /**
     * Functor class representing the objective function of a edge-regularized regression model.
     */
    class edgenet_binomial_likelihood_function
    {
    public:
        /**
         * Creates an objective function object that can be used for minimization using dlib.
         *
         * @param data the complete dataset required for edge-regularized regression
         * @param cvset a cross-validation set
         */
        edgenet_binomial_likelihood_function
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

            arma::Mat<double> B(P, Q);
            // TODO: efficient cast b to B
            for (unsigned int i = 0; i < P; ++i)
                for (unsigned int j = 0; j < Q; ++j)
                    B(i, j) = b(i, j);
            arma::Mat<double> sigm = (X_ * B);
            sigm.transform([](double val) { return log(1 / (1 + exp(val))); });
            // TODO: from here on
            arma::Mat<double> loglik = Y .* log(sigm)  + (1 - Y) .* log(1-sigm);
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
