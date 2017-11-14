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
#include "../../src/graph_functions.hpp"
#include "../../src/edgenet_gaussian.hpp"

class F: public netreg::edgenet_gaussian
{};

/*
* Testing suite for edgenet function for Gaussian variables
*/
BOOST_FIXTURE_TEST_SUITE(edgenet_gaussian_test, F);

BOOST_AUTO_TEST_CASE(test_uucd)
{
    double s = 0;
    double norm = 0;
    const int pi = 0;
    const int qi = 0;
    const double psigx = 0;
    const double psigy = 0;
    const int n = 5;
    const double fill = 1.5;
    const int maxit = 1;

    arma::Mat<double> X(n, n), Y(n, n), B(n, n), GX(n, n), GY(n, n);
    X.fill(fill);
    Y.fill(fill);
    B.fill(fill);
    GX.fill(fill);
    GY.fill(fill);
    for (uint32_t i = 0; i < n; ++i)
    {
      GX(i, i) = 0;
      GY(i, i) = 0;
    }

    arma::Mat<double> lx   = netreg::laplacian(GX.memptr(), n, n, 1);
    arma::Mat<double> ly   = netreg::laplacian(GY.memptr(), n, n, 1);
    arma::rowvec b_row     = B.row(pi);
    arma::rowvec txx_row   = X.row(pi);
    arma::Mat<double> txy  = Y;
    arma::rowvec lx_row    = lx.row(pi);

    set_params(
      s, norm, n, n, pi, qi, psigx, psigy, txx_row, txy,
      lx_row, ly, B, b_row);

    BOOST_REQUIRE(s == netreg::partial_least_squares(txx_row, txy, B, pi, qi));
    BOOST_REQUIRE(norm == txx_row(pi));
}

BOOST_AUTO_TEST_CASE(test_set_param_with_y_penalty)
{
    double s = 1;
    double norm = 1;
    const int pi = 0;
    const int qi = 0;
    const double psigx = 0;
    const double psigy = 1;
    const int n = 5;
    const double fill = 1.5;

    arma::Mat<double> X(n, n), Y(n, n), B(n, n), GX(n, n), GY(n, n);
    X.fill(fill);
    Y.fill(fill);
    B.fill(fill);
    GX.fill(fill);
    GY.fill(fill);
    for (uint32_t i = 0; i < n; ++i)
    {
      GX(i, i) = 0;
      GY(i, i) = 0;
    }

    arma::Mat<double> lx   = netreg::laplacian(GX.memptr(), n, n, 1);
    arma::Mat<double> ly   = netreg::laplacian(GY.memptr(), n, n, 1);
    arma::rowvec b_row     = B.row(pi);
    arma::rowvec txx_row   = X.row(pi);
    arma::Mat<double> txy  = Y;
    arma::rowvec lx_row    = lx.row(pi);

    set_params(
      s, norm, n, n, pi, qi, psigx, psigy, txx_row, txy,
      lx_row, ly, B, b_row);

    BOOST_REQUIRE(s == netreg::partial_least_squares(txx_row, txy, B, pi, qi));
    BOOST_REQUIRE(norm == txx_row(pi));
}

BOOST_AUTO_TEST_CASE(test_set_param_with_x_penalty)
{
    double s = 1;
    double norm = 1;
    const int pi = 0;
    const int qi = 0;
    const double psigx = 1;
    const double psigy = 0;
    const int n = 5;
    const double fill = 1.5;

    arma::Mat<double> X(n, n), Y(n, n), B(n, n), GX(n, n), GY(n, n);
    X.fill(fill);
    Y.fill(fill);
    B.fill(fill);
    GX.fill(fill);
    GY.fill(fill);
    for (uint32_t i = 0; i < n; ++i)
    {
      GX(i, i) = 0;
      GY(i, i) = 0;
    }

    arma::Mat<double> lx   = netreg::laplacian(GX.memptr(), n, n, 1);
    arma::Mat<double> ly   = netreg::laplacian(GY.memptr(), n, n, 1);
    arma::rowvec b_row     = B.row(pi);
    arma::rowvec txx_row   = X.row(pi);
    arma::Mat<double> txy  = Y;
    arma::rowvec lx_row    = lx.row(pi);

    set_params(
      s, norm, n, n, pi, qi, psigx, psigy, txx_row, txy,
      lx_row, ly, B, b_row);

    BOOST_REQUIRE(s == netreg::partial_least_squares(txx_row, txy, B, pi, qi));
    BOOST_REQUIRE(norm == txx_row(pi));
}


BOOST_AUTO_TEST_CASE(test_set_param)
{
    double s = 1;
    double norm = 1;
    const int pi = 0;
    const int qi = 0;
    const double psigx = 0;
    const double psigy = 0;
    const int n = 5;
    const double fill = 1.5;

    arma::Mat<double> X(n, n), Y(n, n), B(n, n), GX(n, n), GY(n, n);
    X.fill(fill);
    Y.fill(fill);
    B.fill(fill);
    GX.fill(fill);
    GY.fill(fill);
    for (uint32_t i = 0; i < n; ++i)
    {
      GX(i, i) = 0;
      GY(i, i) = 0;
    }

    arma::Mat<double> lx   = netreg::laplacian(GX.memptr(), n, n, 1);
    arma::Mat<double> ly   = netreg::laplacian(GY.memptr(), n, n, 1);
    arma::rowvec b_row     = B.row(pi);
    arma::rowvec txx_row   = X.row(pi);
    arma::Mat<double> txy  = Y;
    arma::rowvec lx_row    = lx.row(pi);

    set_params(
      s, norm, n, n, pi, qi, psigx, psigy, txx_row, txy,
      lx_row, ly, B, b_row);

    BOOST_REQUIRE(s == netreg::partial_least_squares(txx_row, txy, B, pi, qi));
    BOOST_REQUIRE(norm == txx_row(pi));
}

BOOST_AUTO_TEST_CASE(test_penalization)
{
    double s = 1;
    double norm = 1;
    const int pi = 0;
    const int qi = 0;
    const int Q = 2;
    const double psigx = 1;
    const double psigy = 1;
    const int m = 5;
    const double fill = 1.5;

    arma::Mat<double> x(m, m);
    arma::Mat<double> cfs(m, m);
    x.fill(fill);
    for (uint32_t i = 0; i < m; ++i) x(i, i) = 0;
    cfs.fill(fill);
    arma::Mat<double> ly = netreg::laplacian(x.memptr(), m, m, 1);
    arma::Mat<double> lx = netreg::laplacian(x.memptr(), m, m, 1);
    arma::rowvec cfs_row = cfs.row(pi);
    arma::rowvec lx_row = lx.row(pi);

    double s_l = s;
    double norm_l = norm;

    lx_penalize(s, norm, pi, qi, psigx, cfs, lx_row);
    ly_penalize(s, norm, pi, qi, psigy, ly, cfs_row);
    graph_penalize(s_l, norm_l, pi, qi, psigx, psigy, Q, lx_row, ly, cfs, cfs_row);
    BOOST_REQUIRE(s_l == s);
    BOOST_REQUIRE(norm_l == norm);
}

BOOST_AUTO_TEST_CASE(test_penalization_only_psigx)
{
    double s = 1;
    double norm = 1;
    const int pi = 0;
    const int qi = 0;
    const int Q = 2;
    const double psigx = 1;
    const double psigy = 0;
    const int m = 5;
    const double fill = 1.5;

    arma::Mat<double> x(m, m);
    arma::Mat<double> cfs(m, m);
    for (uint32_t i = 0; i < m; ++i) x(i, i) = 0;
    x.fill(fill);
    cfs.fill(fill);
    arma::Mat<double> ly = netreg::laplacian(x.memptr(), m, m, 1);
    arma::Mat<double> lx = netreg::laplacian(x.memptr(), m, m, 1);
    arma::rowvec cfs_row = cfs.row(pi);
    arma::rowvec lx_row = lx.row(pi);

    double s_lx = s;
    double norm_lx = norm;

    lx_penalize(s, norm, pi, qi, psigx, cfs, lx_row);
    graph_penalize(s_lx, norm_lx, pi, qi, psigx, psigy, Q, lx_row, ly, cfs, cfs_row);
    BOOST_REQUIRE(s_lx == s);
    BOOST_REQUIRE(norm_lx == norm);
}

BOOST_AUTO_TEST_CASE(test_penalization_only_psigy)
{
    double s = 1;
    double norm = 1;
    const int pi = 0;
    const int qi = 0;
    const int Q = 2;
    const double psigx = 0;
    const double psigy = 1;
    const int m = 5;
    const double fill = 1.5;

    arma::Mat<double> x(m, m);
    arma::Mat<double> cfs(m, m);
    for (uint32_t i = 0; i < m; ++i) x(i, i) = 0;
    x.fill(fill);
    cfs.fill(fill);
    arma::Mat<double> ly = netreg::laplacian(x.memptr(), m, m, 1);
    arma::Mat<double> lx = netreg::laplacian(x.memptr(), m, m, 1);
    arma::rowvec cfs_row = cfs.row(pi);
    arma::rowvec lx_row = lx.row(pi);

    double s_ly = s;
    double norm_ly = norm;

    ly_penalize(s, norm, pi, qi, psigy, ly, cfs_row);
    graph_penalize(s_ly, norm_ly, pi, qi, psigx, psigy, Q,lx_row, ly, cfs, cfs_row);
    BOOST_REQUIRE(s_ly == s);
    BOOST_REQUIRE(norm_ly == norm);
}

BOOST_AUTO_TEST_CASE(test_lx_penalization)
{
    double s = 1;
    double norm = 1;
    int pi = 0;
    int qi = 0;
    const double psigx = 1;
    const int m = 5;
    const double fill = 1.5;
    arma::Mat<double> x(m, m);
    arma::Mat<double> cfs(m, m);
    for (uint32_t i = 0; i < m; ++i) x(i, i) = 0;
    x.fill(fill);
    cfs.fill(fill);
    arma::Mat<double> lx = netreg::laplacian(x.memptr(), m, m, 1);
    arma::rowvec lx_row = lx.row(pi);

    double pen = -lx_row(pi)*cfs(pi, qi) + arma::accu(lx_row  * cfs.col(qi));
    double s_expect = s -  2 * psigx * pen;
    double norm_expect = norm + 2 * psigx * lx_row(pi);

    lx_penalize(s, norm, pi, qi, psigx, cfs, lx_row);
    BOOST_REQUIRE(s_expect == s);
    BOOST_REQUIRE(norm_expect == norm);
}

BOOST_AUTO_TEST_CASE(test_ly_penalization)
{
    double s = 1;
    double norm = 1;
    int pi = 0;
    int qi = 0;
    const double psigy = 1;
    const int m = 5;
    const double fill = 1.5;
    arma::Mat<double> x(m, m);
    arma::rowvec r(m);
    for (uint32_t i = 0; i < m; ++i) x(i, i) = 0;
    x.fill(fill);
    r.fill(fill);
    arma::Mat<double> ly = netreg::laplacian(x.memptr(), m, m, 1);

    double pen = -r(qi) * ly(qi, qi) + arma::accu(r * ly.col(qi));
    double s_expect = s -  2 * psigy * pen;
    double norm_expect = norm + 2 * psigy * ly(qi, qi);

    ly_penalize(s, norm, pi, qi, psigy, ly, r);
    BOOST_REQUIRE(s_expect == s);
    BOOST_REQUIRE(norm_expect == norm);
}

BOOST_AUTO_TEST_SUITE_END()
