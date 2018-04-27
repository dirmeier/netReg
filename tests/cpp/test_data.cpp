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

#include "../../src/graph_functions.hpp"
#include "../../src/stat_functions.hpp"
#include "../../src/family.hpp"
#include "../../src/data_factory.hpp"
#include "../../src/params.hpp"
#include "../../src/edgenet.hpp"
#include "../../src/graph_model_cv_data.hpp"

static double threshold = .0000001;
static int maxit   = 100000;

static double lambda = .1;
static double psi    = .05;
static double phi    = .06;
static uint32_t nfolds = 10;

static uint32_t n = 20;
static uint32_t p = 11;
static uint32_t q = 12;

static std::unique_ptr<double[]> x(new double[n * p]);
static std::unique_ptr<double[]> y(new double[n * q]);
static std::unique_ptr<double[]> gx(new double[p * p]);
static std::unique_ptr<double[]> gy(new double[q * q]);

bool is_identical(arma::Mat<double> &matrix,
                  double *const ptr,
                  uint32_t nrow,
                  uint32_t ncol)
{
    if (matrix.n_rows != nrow || matrix.n_cols != ncol)
    {
        return false;
    }

    for (unsigned int i = 0; i < matrix.n_rows; ++i)
    {
        for (unsigned int j = 0; j < matrix.n_cols; ++j)
        {
            if (matrix(i, j) != ptr[i + ncol * j])
                return false;
        }
    }

    return true;
}

bool is_identical(arma::Mat<double> &m1, arma::Mat<double> &m2)
{
    if (m1.n_rows != m2.n_rows)
        return false;
    if (m2.n_rows != m2.n_cols)
        return false;
    return arma::accu(abs(m1 - m2)) < 0.000001;
}

bool is_identical(arma::rowvec &v1, arma::rowvec &v2)
{
    if (v1.n_elem != v2.n_elem)
        return false;
    return arma::accu(abs(v1 - v2)) < 0.000001;
}

void init_ptrs()
{
    for (uint32_t i = 0; i < n * p; ++i)
        x[i]        = i * 2.0;
    for (uint32_t i = 0; i < n * q; ++i)
        y[i]        = i * 2.0;
    for (uint32_t i = 0; i < p * p; ++i)
        gx[i]       = i * 2.0;
    for (uint32_t i = 0; i < q * q; ++i)
        gy[i]       = i * 2.0;
}

std::map<int, int> count_folds(std::vector<int> &vec)
{
    std::map<int, int> folds;
    for (std::vector<int>::size_type i = 0; i < vec.size(); ++i)
    {
        if (folds.find(vec[i]) == folds.end())
            folds[vec[i]] = 0;
        folds[vec[i]]++;
    }

    return folds;
}


/*
* Testing suite for params
*/
BOOST_AUTO_TEST_SUITE(netReg_parameter_test)

BOOST_AUTO_TEST_CASE(test_params_unset)
{
    netreg::params pars = netreg::params();

    BOOST_REQUIRE(pars.lambda() == 0);
    BOOST_REQUIRE(pars.psigx() == 0);
    BOOST_REQUIRE(pars.psigy() == 0);
    BOOST_REQUIRE(pars.do_lambda() == true);
    BOOST_REQUIRE(pars.do_psigx() == true);
    BOOST_REQUIRE(pars.do_psigy() == true);
    BOOST_REQUIRE(pars.thresh() == .0001);
    BOOST_REQUIRE(pars.niter() == 1000);
    BOOST_REQUIRE(pars.optim_niter() == 1000);
    BOOST_REQUIRE(pars.optim_epsilon() == .0001);
}


BOOST_AUTO_TEST_CASE(test_params_set)
{
    netreg::params pars = netreg::params()
      .lambda(lambda)
      .psigx(psi)
      .psigy(phi)
      .thresh(threshold)
      .niter(maxit)
      .optim_epsilon(threshold)
      .optim_niter(maxit)
      .do_psigx(false)
      .do_psigy(false)
      .do_lambda(false);

    BOOST_REQUIRE(pars.lambda() == lambda);
    BOOST_REQUIRE(pars.psigx() == psi);
    BOOST_REQUIRE(pars.psigy() == phi);
    BOOST_REQUIRE(pars.do_lambda() == false);
    BOOST_REQUIRE(pars.do_psigx() == false);
    BOOST_REQUIRE(pars.do_psigy() == false);
    BOOST_REQUIRE(pars.thresh() == threshold);
    BOOST_REQUIRE(pars.niter() == maxit);
    BOOST_REQUIRE(pars.optim_niter() == maxit);
    BOOST_REQUIRE(pars.optim_epsilon() == threshold);
}

BOOST_AUTO_TEST_SUITE_END()


/*
* Testing suite for the graph cv data class
*/
BOOST_AUTO_TEST_SUITE(netReg_graph_cv_data_test)

BOOST_AUTO_TEST_CASE(test_folds)
{
    init_ptrs();
    std::string fam = "gaussian";
    netreg::graph_model_cv_data dat = netreg::data_factory::build_cv_data(
      x.get(), y.get(), gx.get(), gy.get(), n, p, q, fam, nfolds);

    std::map<int, int> folds = count_folds(dat.fold_ids());
    BOOST_REQUIRE(static_cast<uint32_t>(dat.fold_ids().size()) == n);
    for (auto &entry : folds)
    {
        BOOST_REQUIRE(entry.second == static_cast<int>(n / nfolds));
    }
}

BOOST_AUTO_TEST_CASE(test_cv_set)
{
    init_ptrs();
    std::string fam = "gaussian";
    netreg::graph_model_cv_data dat = netreg::data_factory::build_cv_data(
      x.get(), y.get(), gx.get(), gy.get(), n, p, q, fam, nfolds);

    netreg::cv_set &set    = dat.cvset();
    std::vector<int> folds = dat.fold_ids();

    BOOST_REQUIRE(static_cast<uint32_t>(set.folds().size()) == nfolds);
    for (int i = 0; i < set.fold_count(); ++i)
    {
        netreg::cv_fold &fold = set.get_fold(i);
        for (arma::uvec::iterator j = fold.test_set().begin();
             j != fold.test_set().end(); ++j)
        {
            BOOST_REQUIRE(folds[*j] == i);
        }
        for (arma::uvec::iterator j = fold.train_set().begin();
             j != fold.train_set().end(); ++j)
        {
            BOOST_REQUIRE(folds[*j] != i);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()

/*
* Testing suite for the graph data class
*/
BOOST_AUTO_TEST_SUITE(netReg_graph_data_test)


BOOST_AUTO_TEST_CASE(test_laplacian_matrices)
{
    init_ptrs();
    std::string fam = "gaussian";
    netreg::graph_model_data dat = netreg::data_factory::build_data(
      x.get(), y.get(), gx.get(), gy.get(), n, p, q, fam);

    arma::Mat<double> xm(gx.get(), p, p, false, true);
    arma::Mat<double> ym(gy.get(), q, q, false, true);
    arma::Mat<double> gxx = netreg::laplacian(xm);
    arma::Mat<double> gyy = netreg::laplacian(ym);

    BOOST_REQUIRE(is_identical(dat.lx(), gxx));
    BOOST_REQUIRE(is_identical(dat.ly(), gyy));
}

BOOST_AUTO_TEST_SUITE_END()

/*
* Testing suite for the model data class
*/
BOOST_AUTO_TEST_SUITE(netReg_data_test)

BOOST_AUTO_TEST_CASE(test_dimensions)
{
    init_ptrs();
    arma::Mat<double> xm(x.get(), n, p, false, true);
    arma::Mat<double> ym(y.get(), n, q, false, true);
    netreg::model_data dat = netreg::model_data(
      xm, ym, netreg::family::GAUSSIAN);

    BOOST_REQUIRE(static_cast<uint32_t>(dat.sample_count()) == n);
    BOOST_REQUIRE(static_cast<uint32_t>(dat.response_count()) == q);
    BOOST_REQUIRE(static_cast<uint32_t>(dat.covariable_count()) == p);
}

BOOST_AUTO_TEST_CASE(test_family)
{
    init_ptrs();
    arma::Mat<double> xm(x.get(), n, p, false, true);
    arma::Mat<double> ym(y.get(), n, q, false, true);
    netreg::model_data dat = netreg::model_data(
      xm, ym, netreg::family::GAUSSIAN);


    BOOST_REQUIRE(dat.distribution_family() == netreg::family::GAUSSIAN);
}

BOOST_AUTO_TEST_CASE(test_txx)
{
    init_ptrs();
    arma::Mat<double> xm(x.get(), n, p, false, true);
    arma::Mat<double> ym(y.get(), n, q, false, true);
    netreg::model_data dat = netreg::model_data(
      xm, ym, netreg::family::GAUSSIAN);


    arma::Mat<double> X(x.get(), n, p, false, true);
    arma::Mat<double> txx               = X.t() * X;
    std::vector<arma::rowvec> &txx_rows = dat.txx_rows();

    BOOST_REQUIRE(
      static_cast<uint32_t>(txx.n_rows) == static_cast<uint32_t>(txx_rows.size()));
    for (std::vector<arma::Row<double>>::size_type i = 0; i < txx.n_rows; ++i)
    {
        arma::rowvec v1 = txx.row(i);
        arma::rowvec v2 = txx_rows[i];
        BOOST_REQUIRE(is_identical(v1, v2));
    }
}

BOOST_AUTO_TEST_SUITE_END()
