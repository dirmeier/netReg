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


#ifndef NETREG_CROSSVALIDATION_HPP
#define NETREG_CROSSVALIDATION_HPP

#include <numeric>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

#include "edgenet.hpp"
#include "cv_set.hpp"
#include "graph_model_cv_data.hpp"
#include "error_functions.hpp"

#include "../inst/include/dlib/matrix.h"

namespace netreg
{
    /**
     * Functor class representing a cross validator.
     */
    template<typename T>
    class cross_validator
    {
    public:
        /**
         * Creates an objective function object that can be used for
         * minimization using dlib.
         *
         * @param data the complete dataset required for edge-regularized
         * regression
         *
         * @param cvset a cross-validation set
         */
        cross_validator(
          graph_model_cv_data& data, params& pars):
          pars_(pars),
          cvset_(data.cvset()),
          nfolds_(static_cast<int>(data.cvset().fold_count())),
          model_(data.data(), pars),
          do_lambda_(pars.do_lambda()),
          do_psigx_(pars.do_psigx()),
          do_psigy_(pars.do_psigy())
        {}

        /**
         * Over-write operator () in order to get functor functionality
         *
         * @param params the free parameters on the objective function
         */
        double operator()(const dlib::matrix<double>& p) const
        {
            std::vector<double> sses(nfolds_);
            double lam = do_lambda_  ? p(0, 0): pars_.lambda();
            double psigx = do_psigx_ ? p(1, 0): pars_.psigx();
            double psigy = do_psigy_ ? p(2, 0): pars_.psigy();

            // N-fold cross-validation
            // DO NOT PARALLELIZE.
            // R can't handle it (great)
            for (int fc = 0; fc < nfolds_; ++fc)
            {
                cv_fold& fold = cvset_.get_fold(fc);
                arma::Mat<double> coef = model_.run_cv(
                  lam, psigx, psigy, fold);

                sses[fc] = mse(coef, fold.test_x(), fold.test_y());
            }

            return std::accumulate(sses.begin(), sses.end(), 0.0);
        }

    private:
        params& pars_;
        cv_set& cvset_;     // cv-set on which the selected model is evaluated
        int nfolds_;  // number of folds
        T model_;
        bool do_lambda_;
        bool do_psigx_;
        bool do_psigy_;
    };
}
#endif
