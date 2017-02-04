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

#ifndef NETREG_UTILITY_FUNCTIONS_HPP
#define NETREG_UTILITY_FUNCTIONS_HPP

#include <vector>
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

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
    arma::Mat<double> to_matrix(std::vector< std::vector<double> > &elems);
}

#endif //NETREG_UTILITY_H
