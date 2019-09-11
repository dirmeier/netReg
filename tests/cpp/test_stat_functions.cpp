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

#include "../../src/stat_functions.hpp"

/*
* Testing suite for stats functions
*/
BOOST_AUTO_TEST_SUITE(stat_function_tests)

BOOST_AUTO_TEST_CASE(test_intercept)
{
    int m     = 10;
    double f1 = 2.0;
    double f2 = 1.5;
    double f3 = 3;
    arma::Mat<double> X(m, 1), Y(m, 1), B(1, 1);
    X.fill(f1);
    B.fill(f2);
    Y.fill(f3);
    double expect = f3 - (f2 * f1);
    BOOST_REQUIRE(netreg::intercept(X, Y, B)(0) == expect);
}

BOOST_AUTO_TEST_CASE(test_intercept_non_zero)
{
    int m        = 10;
    double intr  = 3;
    double slope = 2.0;
    arma::Mat<double> X(m, 1), Y(m, 1), B(1, 1);
    B.fill(slope);
    for (int i = 0; i < m; ++i)
    {
        X(i, 0) = i;
        Y(i, 0) = X(i, 0) * B(0, 0) + intr;
    }
    arma::Col<double> intercept = netreg::intercept(X, Y, B);
    BOOST_REQUIRE(intercept(0) == intr);
}

BOOST_AUTO_TEST_CASE(test_intercept_non_zero_matrix)
{
    int m        = 10;
    int n        = 5;
    double intr  = 3;
    double slope = 2.0;
    arma::Mat<double> X(m, n), Y(m, n), B(n, n);
    B.fill(0.0);
    B.diag() += slope;
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            X(i, j) = i;
            Y(i, j) = X(i, j) * B(j, j) + intr;
        }
    }
    arma::Col<double> intercept = netreg::intercept(X, Y, B);
    BOOST_REQUIRE(static_cast<int>(intercept.n_rows) == n);
    BOOST_REQUIRE(intercept(0) == intr);
    BOOST_REQUIRE(intercept(1) == intr);
}

BOOST_AUTO_TEST_CASE(test_pls_zero_indices)
{
    unsigned int p = 2;
    unsigned int q = 2;
    double fill    = 1.0;
    arma::rowvec txx(p);
    for (unsigned int i = 0; i < p; ++i)
    {
        txx(i) = i;
    }
    arma::Mat<double> txy(p, q);
    txy.fill(fill);
    arma::Mat<double> B(p, q);
    B.fill(fill);
    int pi        = 0;
    int qi        = 0;
    double expect = txy(pi, qi) - txx(1) * B(1, qi);
    BOOST_REQUIRE(netreg::partial_least_squares(txx, txy, B, pi, qi) == expect);
}

BOOST_AUTO_TEST_CASE(test_pls_one_indices)
{
    unsigned int p = 2;
    unsigned int q = 2;
    arma::rowvec txx(p);
    arma::Mat<double> txy(p, q);
    arma::Mat<double> B(p, q);
    for (unsigned int i = 0; i < p; ++i)
    {
        txx(i) = i + 2;
        for (unsigned int j = 0; j < q; ++j)
        {
            txy(i, j) = j * i + 1;
            B(i, j)   = -j * i + 3;
        }
    }
    int pi        = 1;
    int qi        = 1;
    double expect = txy(pi, qi) - txx(0) * B(0, qi);
    BOOST_REQUIRE(netreg::partial_least_squares(txx, txy, B, pi, qi) == expect);
}

BOOST_AUTO_TEST_SUITE_END()
