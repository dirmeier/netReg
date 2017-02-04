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

#ifndef NETREG_CV_FOLD_HPP
#define NETREG_CV_FOLD_HPP

#include <vector>
#include <utility>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace netreg
{
/**
 * Class that represents a fold in k-fold cross-validation containing test and training sets.
 */
    class cv_fold
    {
    public:

        cv_fold()
        { }

        cv_fold(std::vector<int> &train_idxs, std::vector<int> &test_idxs):
            train_indexes_(train_idxs.size()),
            test_indexes_(test_idxs.size())
        {
            for (unsigned int j = 0; j < train_idxs.size(); ++j)
                train_indexes_(j) = train_idxs[j];
            for (unsigned int j = 0; j < test_idxs.size(); ++j)
                test_indexes_(j) = test_idxs[j];
        }

        /**
         * Get the indexes of the test set.
         *
         * @return vector of indexes of the test set
         */
        arma::uvec &test_set()
        {
            return test_indexes_;
        }

        /**
         * Get the indexes of the train set.
         *
         * @return vector of indexes of the train set
         */
        arma::uvec &train_set()
        {
            return train_indexes_;
        }

    private:
        arma::uvec train_indexes_; // indexes of train set
        arma::uvec test_indexes_;  // indexes of test set
    };
}
#endif //NETREG_FOLD_HPP
