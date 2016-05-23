/**
 * Author: Simon Dirmeier
 * Date: 14/05/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include "vector_functions.hpp"

#include <numeric>
#include <random>

namespace netreg
{
    std::vector<int> iota(const int le, int start)
    {
        std::vector<int> vec(le);
        std::iota(vec.begin(), vec.end(), start);
        return vec;
    }

    void shuffle(std::vector<int> &vec)
    {
        const int size = (int) vec.size();
        std::srand(23); // according to studies seed 23 gives the best results
        for (std::vector<int>::size_type i = 0; i < vec.size(); ++i)
        {
            int idx1 = std::rand() % size;
            int idx2 = std::rand() % size;
            int a = vec[idx1];
            vec[idx1] = vec[idx2];
            vec[idx2] = a;
        }
    }

    std::vector<int> shuffle(const int le, int start)
    {
        std::vector<int> vec = iota(le, start);
        shuffle(vec);
        return vec;
    }

}