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

#include "cv_set.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <utility>
#include "vector_functions.hpp"

namespace netreg
{
    void cv_set::init(arma::Mat<double> &X, arma::Mat<double> &Y)
    {
        // get a random permutation of the indexes from 1 to N_
        std::vector<int> perms = netreg::shuffle(N_, 0);

        // init fold objects
        std::vector<std::vector<int>> trains(N_FOLDS_);
        std::vector<std::vector<int>> tests(N_FOLDS_);
        for (int i = 0; i < N_FOLDS_; ++i)
        {
            trains.push_back(std::vector<int>());
            tests.push_back(std::vector<int>());
        }

        // running variable
        int test_idx = 0;
        // fill fold objects with sample indexes
        while (test_idx < N_)
        {
            for (int cv = 0; cv < N_FOLDS_; ++cv)
            {
                if (test_idx >= N_)
                    break;
                // add current N_ perms[test_idx] to test set
                if (test_idx < N_)
                {
                    tests[cv].push_back(perms[test_idx]);
                    // add the same index to ALL other train sets
                    for (int k = 0; k < N_FOLDS_; ++k)
                    {
                        if (k != cv)
                            trains[k].push_back(perms[test_idx]);
                    }
                }
                test_idx++;
            }
        }
        //    #pragma omp parallel for
        for (size_t i = 0; i < folds_.size(); ++i)
        {
            folds_[i] = cv_fold(trains[i], tests[i], X, Y);
        }
    }

    /*
     * TODO: this might need debugging ... looks correct though
     */
    void cv_set::init(int *const foldids,
                      arma::Mat<double> &X,
                      arma::Mat<double> &Y)
    {
        throw not_implemented_exception();
    }
    //        for (int i = 0; i < N_FOLDS_; ++i)
    //            folds_.push_back(cv_fold());
    //        // iterate over all sample indeces
    //        for (int i = 0; i < N_; ++i)
    //        {
    //            // get fold id of sample at index i
    //            int fold = foldids[i];
    //            // add the index i to the test set of fold fold
    //            folds_[fold].add_to_test_set(i);
    //            // add index i to all trainings set folds except fold fold
    //            for (int j = 0; j < N_FOLDS_; ++j)
    //            {
    //                if (j != fold) folds_[j].add_to_train_set(i);
    //            }
    //        }
    //    }
}
