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


#ifndef NETREG_OPTIM_HPP
#define NETREG_OPTIM_HPP

#include <string>
#include <vector>
#include <map>
#include <cmath>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#include <iostream>
#endif

#include "graph_model_cv_data.hpp"
#include "cv_set.hpp"
#include "edgenet_gaussian_cv_deviance.hpp"
#include "../inst/include/dlib/optimization.h"

namespace netreg
{
    /**
     * Class that uses several optimization functions from the dlib library.
     * dlib <3
     */
    class optim
    {
    public:
        /**
         * Calculate the pareto optimal points for a given loss function using
         * the BOBYQA algorithm by Powell
         * (MJD Powell, The BOBYQA algorithm for bound constrained optimization
         * without derivatives, 2009)
         *
         * BOBYQA allows for box constraints and does not need derivates.
         * It finds a local optimum for non-convex loss functions defined by a
         * specified radius.
         *
         * @template loss_function the class-name of an objective function that
         *  should be minimized. The class must overwrite the () operator and
         *  return a double.
         * @template deviance the deviance of some distribution
         *
         * @param data the model data for the loss function
         * @param start a vector of starting values for the parameters of the
         *  optimziation
         * @param lower_bound lower bound box constraint
         * @param upper_bound upper bound box constraint
         * @param radius_start initial value of trust region
         * @param radius_stop final value of trust region
         * @param niter maximum calls to the loss function
         *
         */
        template<template<typename ...> class validator, typename deviance>
        std::map<std::string, double> bobyqa(
          graph_model_cv_data& data,
          params& pars,
          std::vector<double>& start,
          std::vector<double>& lower_bound,
          std::vector<double>& upper_bound,
          double radius_start,
          double radius_stop,
          int optim_niter)
        {
            const auto sz = static_cast<int>(start.size());

            // convert to dlib objects
            dlib::matrix<double> par(sz, 1), lb(sz, 1), ub(sz, 1);
            for (int i = 0; i < sz; ++i)
            {
                par(i, 0) = start[i];
                lb(i, 0) = lower_bound[i];
                ub(i, 0) = upper_bound[i];
            }

            // minimize the loss_function
            try
            {
                #ifdef USE_RCPPARMADILLO
                GetRNGstate();
                #endif
                dlib::find_min_bobyqa(loss_function(data, pars),
                                      par,
                                      par.size() * 2 + 1,
                                      lb,
                                      ub,
                                      radius_start,
                                      radius_stop,
                                      optim_niter);
                #ifdef USE_RCPPARMADILLO
                PutRNGstate();
                #endif
            }
            catch (const std::exception& e)
            {
                #ifdef USE_RCPPARMADILLO
                Rprintf("Error estimating optim shrinkage parameters.");
                #else
                std::cerr << "Error estimating optim shrinkage parameters."
                          << std::endl;
                #endif
            }

            return {{"lambda", pars.do_lambda() ? par(0, 0) : pars.lambda()},
                    {"psigx",  pars.do_psigx() ? par(1, 0) : pars.psigx()},
                    {"psigy",  pars.do_psigy() ? par(2, 0) : pars.psigy()}};
        }
    };
}

#endif  // NETREG_OPTIM_HPP
