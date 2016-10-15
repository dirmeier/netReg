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

#ifndef NETREG_EDGENET_BINOMIALLOSS_FUNCTION_HPP
#define NETREG_EDGENET_BINOMIAL_LOSS_FUNCTION_HPP

#include <numeric>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "types.hpp"
#include "edgenet_binomial.hpp"
#include "graph_penalized_linear_model_cv_data.hpp"
#include "cv_set.hpp"
#include "error_functions.hpp"
#include "../inst/include/dlib/matrix.h"

namespace netreg
{
    /**
     * Functor class representing the objective function of a edge-regularized regression model.
     */
    class edgenet_binomial_loss_function
    {
    public:
        /**
         * Creates an objective function object that can be used for minimization using dlib.
         *
         * @param data the complete dataset required for edge-regularized regression
         * @param cvset a cross-validation set
         */
        edgenet_binomial_loss_function
            (graph_penalized_linear_model_cv_data &data):
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
        double operator()(const dlib::matrix<double> &p) const
        {
            std::vector<double> sses(nfolds_);
            // do n-fold cross-validation
//            #pragma omp parallel for
//            for (int fc = 0; fc < nfolds_; ++fc)
//            {
//                cv_fold &fold = cvset_.get_fold(fc);
//                matrix<double> coef;
//                if (do_psigx_ && do_psigy_)
//                {
//                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, p(1, 0),
//                                           p(2, 0), fold);
//                }
//                else if (do_psigy_)
//                {
//                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, 0, p(2, 0),
//                                           fold);
//                }
//                else if (do_psigx_)
//                {
//                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, p(1, 0), 0,
//                                           fold);
//                }
//                else
//                {
//                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, 0, 0, fold);
//                }
//                double err = sse(coef, X_, Y_, fold.test_set());
//                sses[fc] = err;
//            }
//            double  err = std::accumulate(sses.begin(), sses.end(), 0.0);
//            return err;
            return 0.0;
        }

    private:
        // data required for a edge-regularized regression model
        graph_penalized_linear_model_cv_data &data_;
        cv_set &cvset_;          // cv-set on which the selected model is evaluated
        matrix<double> &X_;      // design matrix
        matrix<double> &Y_;      // response matrix
        int nfolds_;             // number of folds
        const edgenet_binomial edgenet_;
        const bool do_psigx_;
        const bool do_psigy_;
    };
}
#endif //NETREG_BINOMIAL_EDGENET_LOSS_FUNCTION_HPP
