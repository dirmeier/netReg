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

#include <cmath>
#include <cstdlib>
#include <iostream>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#include <boost/test/unit_test.hpp>

#include "../../src/error_functions.hpp"


/*
* Testing suite for error functions
*/
BOOST_AUTO_TEST_SUITE(error_function_tests)

BOOST_AUTO_TEST_CASE(test_mse_zero)
{
    int m = 3;
    double fill = 1;
    arma::Mat<double> x(m, m);
    arma::Mat<double> y(m, 1);
    arma::Mat<double> b(m, 1);
    y.fill(2);
    b.fill(1);
    x.fill(0);
    for (int i = 0; i < m; ++i)
    {
      x(i, i) = 2;
    }

    BOOST_REQUIRE(netreg::mse(b, x, y) == 0);
}

BOOST_AUTO_TEST_CASE(test_mse_non_zero)
{
    int m = 3;
    arma::Mat<double> x(m, m);
    arma::Mat<double> y(m, 1);
    arma::Mat<double> b(m, 1);
    y.fill(2);
    b.fill(1);
    x.fill(0);
    for (int i = 0; i < m; ++i)
    {
      x(i, i) = 1;
    }

    double expect = ((2 - 1) * m ) * ((2-1) * m) / m;
    BOOST_REQUIRE(netreg::mse(b, x, y) == expect);
}

BOOST_AUTO_TEST_SUITE_END()
