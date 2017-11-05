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

#include "../../src/graph_functions.hpp"

bool degrees_are_correct(std::vector<double>& degrees, arma::Mat<double>& x)
{
    for (std::vector<double>::size_type i = 0; i < degrees.size(); ++i)
    {
        double rowsum = 0;
        for (uint32_t j = 0; j < x.n_cols; ++j)
            rowsum += x(i, j);
        if (degrees[i] != rowsum)
            return false;
    }

    return true;
}

bool laplacian_is_correct(arma::Mat<double>& lapl, arma::Mat<double>& x)
{
    std::vector<double> degrees =
      netreg::degree_distribution(x.memptr(), x.n_rows, x.n_cols);
    for (uint32_t i = 0; i < lapl.n_rows; ++i)
    {
        for (uint32_t j = 0; j < lapl.n_cols; ++j)
        {
            if (lapl(i, j) != 1 - (x(i, j) / degrees[i]) &&
                lapl(i, j) != -x(i, j) / std::sqrt(degrees[i] * degrees[j]) &&
                lapl(i, j) != 0)
            {
                return false;
            }
        }
    }
    return true;
}

/*
* Testing suite for graph functions
*/
BOOST_AUTO_TEST_SUITE(graph_function_tests)

BOOST_AUTO_TEST_CASE(test_degree_distribution)
{
    int m = 3;
    arma::Mat<double> x(m, m);
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            if (i == j)
                continue;
            x(i, j) = i + j;
        }
    }
    std::vector<double> degrees = netreg::degree_distribution(x.memptr(), m, m);
    BOOST_REQUIRE(degrees_are_correct(degrees, x));
}

BOOST_AUTO_TEST_CASE(test_laplacian)
{
    int m = 3;
    arma::Mat<double> x(m, m);
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            if (i == j)
                x(i, j) = 0;
            else
                x(i, j) = i + j;
        }
    }
    arma::Mat<double> lapl = netreg::laplacian(x.memptr(), m, m, 1);
    BOOST_REQUIRE(laplacian_is_correct(lapl, x));
}

BOOST_AUTO_TEST_SUITE_END()
