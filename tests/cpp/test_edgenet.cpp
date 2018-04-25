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

#include "../../src/params.hpp"
#include "../../src/stat_functions.hpp"
#include "../../src/graph_functions.hpp"
#include "../../src/graph_model_data.hpp"
#include "../../src/data_factory.hpp"
#include "../../src/edgenet.hpp"


/*
* Testing suite for edgenet function for Gaussian variables
*/
BOOST_AUTO_TEST_SUITE(edgenet_test);

    BOOST_AUTO_TEST_CASE(test_set_param_with_y_penalty)
    {
        double s = 1;
        double no = 1;
        const int pi = 0;
        const int qi = 0;
        const double lambda = 3;
        const double psigx = 0;
        const double psigy = 1;
        const int n = 5;
        double fill = 1;
        int niter = 100;
        double thresh = 100;

        std::string fam = "gaussian";
        arma::Mat<double> X(n, n);
        arma::Mat<double> Y(n, n);
        arma::Mat<double> B(n, n);
        arma::Mat<double> GX(n, n);
        arma::Mat<double> GY(n, n);
        X.fill(fill);
        Y.fill(fill);
        B.fill(fill);
        GX.fill(fill);
        GY.fill(fill);

        for (int i = 0; i < n; ++i)
        {
            GX(i, i) = 0;
            GY(i, i) = 0;
        }

        arma::Mat<double> lx = netreg::laplacian(GX);
        arma::Mat<double> ly = netreg::laplacian(GY);
        arma::rowvec b_row = B.row(pi);
        arma::rowvec txx_row = X.row(pi);
        arma::Mat<double> txy = Y;
        arma::rowvec lx_row = lx.row(pi);

        netreg::params pars = netreg::params()
          .lambda(lambda)
          .psigx(psigx)
          .psigy(psigy)
          .thresh(thresh)
          .niter(niter);
        netreg::graph_model_data dat = netreg::data_factory::build_data(
          X.memptr(), Y.memptr(), GX.memptr(), GY.memptr(), n, n, n, fam);
        netreg::edgenet edge(dat, pars);

        s = edge.partial(pi, qi, txx_row, txy, B, b_row);
        no = edge.norm(pi, qi, txx_row);

        BOOST_REQUIRE(
          s == netreg::partial_least_squares(txx_row, txy, B, pi, qi)
               -2 * psigy *
               (-b_row(qi) * ly(qi, qi) +
                arma::accu(b_row * ly.col(qi))));
        BOOST_REQUIRE(no == txx_row(pi) + 2 * psigy * ly(qi, qi));
    }

    BOOST_AUTO_TEST_CASE(test_set_param_with_x_penalty)
    {
        double s = 1;
        double no = 1;
        const int pi = 0;
        const int qi = 0;
        const double lambda = 3;
        const double psigx = 1;
        const double psigy = 0;
        const int n = 5;
        double fill = 1;
        int niter = 100;
        double thresh = 100;

        std::string fam = "gaussian";
        arma::Mat<double> X(n, n);
        arma::Mat<double> Y(n, n);
        arma::Mat<double> B(n, n);
        arma::Mat<double> GX(n, n);
        arma::Mat<double> GY(n, n);
        X.fill(fill);
        Y.fill(fill);
        B.fill(fill);
        GX.fill(fill);
        GY.fill(fill);

        for (int i = 0; i < n; ++i)
        {
            GX(i, i) = 0;
            GY(i, i) = 0;
        }

        arma::Mat<double> lx = netreg::laplacian(GX);
        arma::Mat<double> ly = netreg::laplacian(GY);
        arma::rowvec b_row = B.row(pi);
        arma::rowvec txx_row = X.row(pi);
        arma::Mat<double> txy = Y;
        arma::rowvec lx_row = lx.row(pi);

        netreg::params pars = netreg::params()
          .lambda(lambda)
          .psigx(psigx)
          .psigy(psigy)
          .thresh(thresh)
          .niter(niter);
        netreg::graph_model_data dat = netreg::data_factory::build_data(
          X.memptr(), Y.memptr(), GX.memptr(), GY.memptr(), n, n, n, fam);
        netreg::edgenet edge(dat, pars);

        s = edge.partial(pi, qi, txx_row, txy, B, b_row);
        no = edge.norm(pi, qi, txx_row);

        BOOST_REQUIRE(
          s == netreg::partial_least_squares(txx_row, txy, B, pi, qi)
                 -2 * psigx *
                 (-lx(pi,pi) * B(pi, qi) +
                  arma::accu(lx_row * B.col(qi))));
        BOOST_REQUIRE(no == txx_row(pi) + 2 * psigx * lx(pi,pi));
    }

    BOOST_AUTO_TEST_CASE(test_set_param)
    {
        double s = 1;
        double no = 1;
        const int pi = 0;
        const int qi = 0;
        const double lambda = 3;
        const double psigx = 1;
        const double psigy = 0;
        const int n = 5;
        double fill = 1;
        int niter = 100;
        double thresh = 100;

        std::string fam = "gaussian";
        arma::Mat<double> X(n, n);
        arma::Mat<double> Y(n, n);
        arma::Mat<double> B(n, n);
        arma::Mat<double> GX(n, n);
        arma::Mat<double> GY(n, n);
        X.fill(fill);
        Y.fill(fill);
        B.fill(fill);
        GX.fill(fill);
        GY.fill(fill);

        for (int i = 0; i < n; ++i)
        {
            GX(i, i) = 0;
            GY(i, i) = 0;
        }

        arma::Mat<double> lx = netreg::laplacian(GX);
        arma::Mat<double> ly = netreg::laplacian(GY);
        arma::rowvec b_row = B.row(pi);
        arma::rowvec txx_row = X.row(pi);
        arma::Mat<double> txy = Y;
        arma::rowvec lx_row = lx.row(pi);

        netreg::params pars = netreg::params()
          .lambda(lambda)
          .psigx(psigx)
          .psigy(psigy)
          .thresh(thresh)
          .niter(niter);
        netreg::graph_model_data dat = netreg::data_factory::build_data(
          X.memptr(), Y.memptr(), GX.memptr(), GY.memptr(), n, n, n, fam);
        netreg::edgenet edge(dat, pars);

        s = edge.partial(pi, qi, txx_row, txy, B, b_row);
        no = edge.norm(pi, qi, txx_row);

        BOOST_REQUIRE(
          s == netreg::partial_least_squares(txx_row, txy, B, pi, qi));
        BOOST_REQUIRE(no == txx_row(pi));
    }


BOOST_AUTO_TEST_SUITE_END()
