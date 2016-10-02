/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_CVSET_HPP
#define NETREG_CVSET_HPP

#include <vector>

#include "cv_fold.hpp"
#include "math_functions.hpp"

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
        cv_set(const int n, const int n_folds):
            N_FOLDS_(n_folds), N_(n)
        {
            init();

        }

        /**
         * Creates a cross-validation set.
         *
         * @param size the number of samples provided (e.g. 100)
         * @param foldids fold assignments for all samples
         */
        cv_set(const int n, int *const foldids): N_FOLDS_(n), N_(n)
        {
            init(foldids);
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
         * @return the numer of samples
         */
        const int n()
        {
            return N_;
        }

    private:
        // init folds from scratch
        void init();
        // init folds using predefined fold ids
        void init(int *const foldids);

        const int N_FOLDS_;           // the number of folds
        const int N_;              // the sample size
        std::vector <cv_fold> folds_; // the fold objects
    };
}
#endif //NETREG_CVSET_HPP
