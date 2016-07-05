/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "cv_set.hpp"
#include "vector_functions.hpp"

void cv_set::init()
{
    // get a random permutation of the indexes from 1 to SIZE
    std::vector<int> perms = netreg::shuffle(SIZE_, 0);
    // init fold objects
    for (int i = 0; i < N_FOLDS_; ++i)
        folds_.push_back(cv_fold());
    // running variable
    int test_idx = 0;
    // fill fold objects with sample indexes
    while (test_idx < SIZE_)
    {
        for (int cv = 0; cv < N_FOLDS_; ++cv)
        {
            if (test_idx >= SIZE_)
                break;
            // add current index perms[test_idx] to test set
            if (test_idx < SIZE_)
                folds_[cv].add_to_test_set(perms[test_idx]);
            // add the same index to ALL other train sets
            for (int k = 0; k < N_FOLDS_; ++k)
            {
                if (k == cv)
                    continue;
                if (test_idx < SIZE_)
                    folds_[k].add_to_train_set(perms[test_idx]);
            }
            test_idx++;
        }
    }
}
