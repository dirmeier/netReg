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
#include "../../src/family.hpp"
#include "../../src/edgenet_gaussian.hpp"
#include "../../src/graph_penalized_linear_model_cv_data.hpp"
#include "../../src/edgenet_gaussian_model_selection.hpp"

static double threshold = .0000001;
static uint32_t maxit = 100000;

static double lambda = .1;
static double psi = .0;
static double phi = .0;

static uint32_t nfolds = 10;

static uint32_t n = 20;
static uint32_t p = 11;
static uint32_t q = 12;

static std::unique_ptr<double[]> x(new double[n * p]);
static std::unique_ptr<double[]> y(new double[n * q]);
static std::unique_ptr<double[]> gx(new double[p * p]);
static std::unique_ptr<double[]> gy(new double[q * q]);

void init_ptrs()
{
  for (uint32_t i = 0; i < n * p; ++i)
    x[i] = 1;
  for (uint32_t i = 0; i < n * q; ++i)
    y[i] = 1;
  for (uint32_t i = 0; i < p * p; ++i)
    gx[i] = 1;
  for (uint32_t i = 0; i < q * q; ++i)
    gy[i] = 1;
}

std::map<int, int> count_folds(std::vector<int>& vec)
{
  std::map<int, int> folds;
  for (std::vector<int>::size_type i = 0; i < vec.size(); ++i)
  {
    if (folds.find(vec[i]) == folds.end()) folds[vec[i]] = 0;
    folds[vec[i]]++;
  }

  return folds;
}

/*
* Testing suite for the graph_penalized_linear_model_cv_data class
*/
BOOST_AUTO_TEST_SUITE(netReg_cv_data_test)

BOOST_AUTO_TEST_CASE(test_folds)
{
    init_ptrs();
    netreg::graph_penalized_linear_model_cv_data dat = netreg::graph_penalized_linear_model_cv_data(
      x.get(), y.get(), gx.get(), gy.get(),
      n, p, q, lambda, 0, psi, phi, maxit, threshold, nfolds,
      netreg::family::GAUSSIAN);

    std::map<int, int> folds = count_folds(dat.fold_ids());
    BOOST_REQUIRE(static_cast<uint32_t>(dat.fold_ids().size()) == n);
    for (auto& entry: folds)
    {
        BOOST_REQUIRE(entry.second == static_cast<int>(n / nfolds));
    }
}

BOOST_AUTO_TEST_CASE(test_cv_set)
{
    init_ptrs();
    netreg::graph_penalized_linear_model_cv_data dat = netreg::graph_penalized_linear_model_cv_data(
      x.get(), y.get(), gx.get(), gy.get(),
      n, p, q, lambda, 0, psi, phi, maxit, threshold, nfolds,
      netreg::family::GAUSSIAN);

    netreg::cv_set& set = dat.cvset();
    std::vector<int> folds = dat.fold_ids();

    BOOST_REQUIRE(static_cast<uint32_t>(set.folds().size()) == nfolds);
    for (int i = 0; i < set.fold_count(); ++i)
    {
        netreg::cv_fold& fold = set.get_fold(i);
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
* Testing suite for the graph_penalized_linear_model_data class
*/
BOOST_AUTO_TEST_SUITE(netReg_data_test)

BOOST_AUTO_TEST_SUITE_END()
