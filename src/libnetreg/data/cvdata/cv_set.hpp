#ifndef NETREG_CVSET_HPP
#define NETREG_CVSET_HPP

#include <vector>
#include "cv_fold.hpp"

/**
 * Class that respresents a k-fold cross-validation set.
 * Thus k folds are included in one cv_set object.
 */
class cv_set
{
public:

    /**
     * Creates a cross-validation set.
     *
     * @param size the number of samples provided (e.g. 100)
     * @param folds the number of folds to be created (e.g. 10)
     */
    cv_set(const int size, const int folds) :
        N_FOLDS_(folds), SIZE_(size)
    {
        init();
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
    const int size()
    {
        return SIZE_;
    }

private:
    /**
     * Initialize the cross-validation folds
     */
    void init();

    const int N_FOLDS_;           // the number of folds
    const int SIZE_;              // the sample size
    std::vector <cv_fold> folds_; // the fold objects
};

#endif //NETREG_CVSET_HPP
