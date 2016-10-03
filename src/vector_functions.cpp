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
        // TODO: change to R random
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