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

#ifndef NETREG_VECTOR_FUNCTIONS_HPP
#define NETREG_VECTOR_FUNCTIONS_HPP

#include <vector>

namespace netreg
{
    /**
     * Generate a vector of length le with increasing values starting from
     * start.
     *
     * @param le length of vector
     * @param start starting value at index 0 of the newly created vector
     *
     * @return a vector with increasing values
     */
    std::vector<int> iota(int le, int start);

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
    std::vector<int> shuffle(int le, int start);
}
#endif  // NETREG_VECTOR_FUNCTIONS_HPP
