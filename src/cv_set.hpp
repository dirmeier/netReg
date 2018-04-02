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

#ifndef NETREG_CV_SET_HPP
#define NETREG_CV_SET_HPP

#include <vector>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif

#include "cv_fold.hpp"
#include "math_functions.hpp"

#include "not_implemented_exception.hpp"

/**
 * Class that respresents a k-fold cross-validation set.
 * Thus k folds are included in one cv_set object.
 */
namespace netreg
{
    class cv_set
    {
       public:
        /**
         * Creates a cross-validation set.
         *
         * @param n the number of samples provided (e.g. 100)
         * @param n_folds the number of folds to be created (e.g. 10)
         */
        cv_set(const int n,
               const int n_folds,
               arma::Mat<double> &X,
               arma::Mat<double> &Y)
            : N_FOLDS_(n_folds), N_(n), folds_(N_FOLDS_)
        {
            init(X, Y);
        }

        /**
         * Creates a cross-validation set.
         *
         * @param size the number of samples provided (e.g. 100)
         * @param foldids fold assignments for all samples
         */
        cv_set(const int n,
               int *const foldids,
               arma::Mat<double> &X,
               arma::Mat<double> &Y)
            : N_FOLDS_(n), N_(n)
        {
            throw not_implemented_exception();
            //init(foldids, X, Y);
        }

        std::vector<cv_fold> &folds()
        {
            return folds_;
        }

        /**
         * Getter for the i-th fold.
         *
         * @param i index of fold
         * @return a cv_fold object
         */
        cv_fold &get_fold(const int i)
        {
            return folds_[i];
        }

        /**
         * Getter for the number of folds.
         *
         * @return the number of folds
         */
        const int fold_count()
        {
            return N_FOLDS_;
        }

        /**
         * Getter for the sample size
         *
         * @return the number of samples
         */
        const int n()
        {
            return N_;
        }

       private:
        // init folds from scratch
        void init(arma::Mat<double> &X, arma::Mat<double> &Y);
        // init folds using predefined fold ids
        void init(int *const foldids,
                  arma::Mat<double> &X,
                  arma::Mat<double> &Y);

        const int N_FOLDS_;           // the number of folds
        const int N_;                 // the sample size
        std::vector<cv_fold> folds_;  // the fold objects
    };
}
#endif  // NETREG_CVSET_HPP
