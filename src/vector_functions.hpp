/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_VECTOR_FUNCTIONS_HPP
#define NETREG_VECTOR_FUNCTIONS_HPP

#include <vector>

namespace netreg
{
    /**
     * Generate a vector of length le with increasing values starting from start.
     *
     * @param le length of vector
     * @param start starting value at index 0
     * @return a vector with increasing values
     */
    std::vector<int> iota(const int le, int start);

    /**
     * Randomize the values of a vector.
     *
     * @param vec a vector of ints
     */
    void shuffle(std::vector<int> &vec);

    /**
     * Generate a vector of length le with increasing values starting from start
     * in random order
     *
     * @param le length of vector
     * @param start starting value at index 0
     * @return a vector with increasing values in random order
     */
    std::vector<int> shuffle(const int le, int start);
}
#endif //NETREG_VECTOR_FUNCTIONS_HPP
