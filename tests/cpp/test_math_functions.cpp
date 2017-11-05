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

#include "../../src/math_functions.hpp"

/*
* Testing suite for math functions
*/
BOOST_AUTO_TEST_SUITE(math_function_tests)

BOOST_AUTO_TEST_CASE(test_abs_dprod)
{
    uint32_t n   = 100;
    double afill = 1;
    double bfill = 2;
    arma::Col<double> a(n);
    a.fill(afill);
    arma::Col<double> b(n);
    b.fill(bfill);
    double pd = netreg::abs_dprod(a, b);
    BOOST_REQUIRE(pd == afill * bfill * n);
}

BOOST_AUTO_TEST_CASE(test_abs_dprod_negative)
{
    uint32_t n   = 100;
    double afill = 1;
    double bfill = -2;
    arma::Col<double> a(n);
    a.fill(afill);
    arma::Col<double> b(n);
    b.fill(bfill);
    double pd = netreg::abs_dprod(a, b);
    BOOST_REQUIRE(pd == (-1) * afill * bfill * n);
}

BOOST_AUTO_TEST_CASE(test_max_element_int)
{
    int n    = 10;
    int *ptr = new int[n];
    for (int i = 0; i < n; ++i)
        ptr[i] = i;
    BOOST_REQUIRE(netreg::max_element(ptr, n) == n - 1);
    delete[] ptr;
}

BOOST_AUTO_TEST_CASE(test_max_element_double)
{
    int n       = 10;
    double *ptr = new double[n];
    for (int i = 0; i < n; ++i)
        ptr[i] = i;
    BOOST_REQUIRE(netreg::max_element(ptr, n) == n - 1);
    delete[] ptr;
}

BOOST_AUTO_TEST_CASE(test_softnorm_full_penalty)
{
    double s     = 1.0;
    double lalph = 1.0;
    double norm  = 1.0;
    BOOST_REQUIRE(netreg::softnorm(s, lalph, norm) == 0);
}

BOOST_AUTO_TEST_CASE(test_softnorm_lasso_pos_penalty)
{
    double s     = 2.0;
    double lalph = 1.0;
    double norm  = 1.0;
    BOOST_REQUIRE(netreg::softnorm(s, lalph, norm) == 1.0);
}

BOOST_AUTO_TEST_CASE(test_softnorm_lasso_neg_penalty)
{
    double s     = -2.0;
    double lalph = 1.0;
    double norm  = 1.0;
    BOOST_REQUIRE(netreg::softnorm(s, lalph, norm) == -1.0);
}

BOOST_AUTO_TEST_CASE(test_softnorm_lasso_pos_penalty_twice)
{
    double s     = 3.0;
    double lalph = 2.0;
    double norm  = 2.0;
    BOOST_REQUIRE(netreg::softnorm(s, lalph, norm) == .5);
}

BOOST_AUTO_TEST_SUITE_END()
