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

#include <cstdlib>
#include <map>
#include <iostream>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#include <boost/test/unit_test.hpp>

#include "../../src/vector_functions.hpp"

/*
 * Testing suite for vector functions
 */
BOOST_AUTO_TEST_SUITE(vector_function_tests)

BOOST_AUTO_TEST_CASE(test_iota)
{
    std::vector<int> v   = netreg::iota(10, 0);
    int              run = 0;
    for (int i : v)
    {
        BOOST_REQUIRE(i == run++);
    }
}

BOOST_AUTO_TEST_CASE(test_shuffle)
{
    std::vector<int> v = netreg::iota(10, 0);
    netreg::shuffle(v);
    int has_shuffle = 0;
    for (std::vector<int>::size_type i = 0; i < v.size() - 1; ++i)
    {
        has_shuffle += v[i] > v[i + 1] ? 1 : 0;
    }
    BOOST_REQUIRE(has_shuffle > 0);
}

BOOST_AUTO_TEST_SUITE_END()
