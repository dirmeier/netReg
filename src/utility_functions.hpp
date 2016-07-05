/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_UTILITY_H
#define NETREG_UTILITY_H

#include <vector>
#include <string>

#include "types.hpp"

namespace netreg
{
    /**
     * Split a string at a specific char and parse to a double vectpr
     *
     * @param str the string you wanna split
     * @param spl the split character
     * @return returns a vector of parsed double values
     */
    std::vector<double> to_double(char* str, const char * spl=",");

    /**
    * Split a string at a specific char and parse to a double vectpr
    *
    * @param str the string you wanna split
    * @param delim the split string
    * @return returns a vector of parsed double values
    */
    std::vector<double> to_double(const std::string& str, const char delim);

    /**
     * Reads a vector of vector of doubles and converts it to a matrix
     *
     * @param elems a vector of a vector of doubles
     * @return a pointer to double
     */
    matrix<double> to_matrix(std::vector< std::vector<double> > &elems);
}

#endif //NETREG_UTILITY_H
