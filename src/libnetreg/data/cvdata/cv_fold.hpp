#ifndef NETREG_FOLD_HPP
#define NETREG_FOLD_HPP

#include <vector>

/**
 * Class that represents a fold in k-fold cross-validation containing test and training sets.
 */
class cv_fold
{
public:
    /**
     * Add an index/sample to the test set
     *
     * @param i the index of the sample
     */
    void add_to_test_set(const int i)
    {
        test_indexes_.push_back(i);
    }

    /**
     * Add an index/sample to the train set
     *
     * @param i the index of the sample
     */
    void add_to_train_set(const int i)
    {
        train_indexes_.push_back(i);
    }

    /**
     * Get the indexes of the test set.
     *
     * @return vector of indexes of the test set
     */
    std::vector<int> &test_set()
    {
        return test_indexes_;
    }

    /**
     * Get the indexes of the train set.
     *
     * @return vector of indexes of the train set
     */
    std::vector<int> &train_set()
    {
        return train_indexes_;
    }

private:
    std::vector<int> train_indexes_; // indexes of train set
    std::vector<int> test_indexes_;  // indexes of test set
};

#endif //NETREG_FOLD_HPP
