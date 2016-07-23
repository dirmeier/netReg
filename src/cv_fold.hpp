/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_FOLD_HPP
#define NETREG_FOLD_HPP

#include <vector>
#include <utility>
#include <iostream>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

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
            for (int j = 0; j < train_idxs.size(); ++j)
                train_indexes_(j) = train_idxs[j];
            for (int j = 0; j < test_idxs.size(); ++j)
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
