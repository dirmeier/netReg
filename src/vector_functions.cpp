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

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include <random>
#include "armadillo"
#endif

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
        #ifdef USE_RCPPARMADILLO
        GetRNGstate();
        Rcpp::Environment base_env("package:base");
        Rcpp::Function set_seed_r = base_env["set.seed"];
        set_seed_r(23);
        #endif
        for (std::vector<int>::size_type i = 0; i < vec.size(); ++i)
        {
            #ifdef USE_RCPPARMADILLO
            int idx = unif_rand() * size;
            #else
            int idx = std::rand() % size;
            #endif
            int a = vec[idx];
            vec[idx] = vec[i];
            vec[i] = a;
        }
        #ifdef USE_RCPPARMADILLO
        PutRNGstate();
        #endif
    }

    std::vector<int> shuffle(const int le, int start)
    {
        std::vector<int> vec = iota(le, start);
        shuffle(vec);
        return vec;
    }

}
