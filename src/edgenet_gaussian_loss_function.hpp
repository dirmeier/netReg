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

#ifndef NETREG_EDGENET_GAUSSIAN_LOSSFUNCTION_HPP
#define NETREG_EDGENET_GAUSSIAN_LOSSFUNCTION_HPP

#include <numeric>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#ifndef NDEBUG
#include <iostream>
#endif
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "edgenet_gaussian.hpp"
#include "graph_penalized_linear_model_cv_data.hpp"
#include "cv_set.hpp"
#include "error_functions.hpp"

#include "../inst/include/dlib/matrix.h"

namespace netreg
{
    /**
     * Functor class representing the objective function of a edge-regularized regression model.
     */
    class edgenet_gaussian_loss_function
    {
    public:
        /**
         * Creates an objective function object that can be used for minimization using dlib.
         *
         * @param data the complete dataset required for edge-regularized regression
         * @param cvset a cross-validation set
         */
        edgenet_gaussian_loss_function
            (graph_penalized_linear_model_cv_data& data):
            data_(data),
            cvset_(data.cvset()),
            X_(data.design()),
            Y_(data.response()),
            nfolds_(static_cast<int>(data.cvset().fold_count())),
            edgenet_(),
            do_psigx_(data.psigx() == -1),
            do_psigy_(data.psigy() == -1)
        { }

        /**
         * Over-write operator () in order to get functor functionality (object behaves like a function)
         *
         * @param params the free parameters on the objective function
         */
        double operator()(const dlib::matrix<double>& p) const
        {
            std::vector<double> sses(nfolds_);
            // do n-fold cross-validation
            #pragma omp parallel for
            for (int fc = 0; fc < nfolds_; ++fc)
            {
                cv_fold& fold = cvset_.get_fold(fc);
                arma::Mat<double> coef;
                if (do_psigx_ && do_psigy_)
                {
                    #ifndef NDEBUG
                    std::cout << "Doing lambda, psigx, psigy: " <<  p(0, 0) << ", " <<  p(1, 0) <<  ", "  << p(2, 0) << std::endl:
                    #endif
                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, p(1, 0), p(2, 0), fold);
                }
                else if (do_psigy_)
                {
                    #ifndef NDEBUG
                    std::cout << "Doing lambda, psigy: " <<  p(0, 0) <<  " "  << p(2, 0) << std::endl:
                    #endif
                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, 0, p(2, 0), fold);
                }
                else if (do_psigx_)
                {
                    #ifndef NDEBUG
                    std::cout << "Doing lambda, psigx: " <<  p(0, 0) <<  " "  << p(1, 0) << std::endl:
                    #endif
                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, p(1, 0), 0, fold);
                }
                else
                {
                    #ifndef NDEBUG
                    std::cout << "Doing lambda: " <<  p(0, 0) << std::endl:
                    #endif
                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, 0, 0, fold);
                }
                double err = sse(coef, fold.test_x(), fold.test_y());
                sses[fc] = err;
            }
            double err = std::accumulate(sses.begin(), sses.end(), 0.0);
            return err;
        }

    private:
        // data required for a edge-regularized regression model
        graph_penalized_linear_model_cv_data& data_;
        cv_set& cvset_;          // cv-set on which the selected model is evaluated
        arma::Mat<double>& X_;      // design matrix
        arma::Mat<double>& Y_;      // response matrix
        const int nfolds_;             // number of folds
        const edgenet_gaussian edgenet_;
        const bool do_psigx_;
        const bool do_psigy_;
    };
}
#endif //NETREG_EDGENETLOSSFUNCTION_HPP
