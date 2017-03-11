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

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "graph_penalized_linear_model_cv_data.hpp"
#include "cv_set.hpp"
#include "edgenet_gaussian_loss_function.hpp"

#include "graph_penalized_linear_model_cv_data.hpp"

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
         * Calculate the pareto optimal points for a given loss function using the BOBYQA algorithm by Powell
         * (MJD Powell, The BOBYQA algorithm for bound constrained optimization without derivatives, 2009)
         *
         * BOBYQA allows for box constraints and does not need derivates.
         * It finds a local optimum for non-convex loss functions defined by a specified radius.
         *
         * We use it to minimize functions from the loss_functions package
         *
         * @template loss_function the class-name of an objective function that should be minimized
         * @param point an intial parameter setting for bobyqa
         * @param data the model data for the loss function
         * @param cvset the cv-set on which the estimated predictor is tested against
         * @param lower_bound lower bound box constraint
         * @param upper_bound upper bound box constraint
         * @param radius_start initial value of trust region
         * @param radius_start final value of trust region
         * @param niter maximum calls to the loss function
         */
        template<typename loss_function>
        std::map<std::string, double> bobyqa
            (graph_penalized_linear_model_cv_data &data,
             std::vector<double> &start,
             std::vector<double> &lower_bound,
             std::vector<double> &upper_bound,
             const double radius_start,
             const double radius_stop,
             const int niter) const
             {
                 const int sz = static_cast<int>(start.size());
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
                     dlib::find_min_bobyqa(
                         loss_function(data),
                         par,
                         par.size() * 2 + 1,
                         lb,
                         ub,
                         radius_start,
                         radius_stop,
                         niter);
                     #ifdef USE_RCPPARMADILLO
                     PutRNGstate();
                     #endif
                 }
                 catch (const std::exception &e)
                 {
                     #ifdef USE_RCPPARMADILLO
                     Rprintf("Error estimating optim shrinkage parameters.");
                     #else
                     std::cerr << "Error estimating optim shrinkage parameters." << std::endl;
                     #endif
                 }
                 return {{"lambda", par(0, 0)},
                         {"psigx",  data.psigx() == -1 ? par(1, 0) : 0.0},
                         {"psigy",  data.psigy() == -1 ? par(2, 0) : 0.0}};
             }

        // todo
        template<typename loss_function>
        std::map<std::string, double> bifurcation
          (graph_penalized_linear_model_cv_data &data,
           std::vector<double> &lower_bound,
           std::vector<double> &upper_bound,
           const double epsilon,
           const int niter) const
           {
             const int sz = static_cast<int>(lower_bound.size());
             loss_function loss(data);
             // convert to dlib objects
             dlib::matrix<double> par(sz, 1);
             for (int i = 0; i < sz; ++i) par(i, 0) = lower_bound[i];
             // minimize the loss_function
             try
             {
                for (std::vector<double>::size_type i = 0;
                     i < upper_bound.size();
                     ++i)
                  {
                      double opt = blockwise_bifurcation<loss_function>(
                        par, i, loss,
                        lower_bound, upper_bound,
                        epsilon, niter);
                      par(i, 0) = opt;
                      break;
                  }
             }
             catch (const std::exception &e)
             {
                 #ifdef USE_RCPPARMADILLO
                 Rprintf("Error estimating optim shrinkage parameters.");
                 #else
                 std::cerr << "Error estimating optim shrinkage parameters." << std::endl;
                 #endif
             }
             return {{"lambda", par(0, 0)},
                     {"psigx",  data.psigx() == -1 ? par(1, 0) : 0.0},
                     {"psigy",  data.psigy() == -1 ? par(2, 0) : 0.0}};
           }

    private:

        template<typename loss_function>
        double blockwise_bifurcation
          (dlib::matrix<double>& par, const int idx,
           loss_function &loss,
           std::vector<double> &lower_bound,
           std::vector<double> &upper_bound,
           const double epsilon,
           const int niter) const
           {

               int iter = 0;
               double l_left  = lower_bound[idx];
               double l_right = upper_bound[idx];
               double err_old = 100000.0;
               double err_new = 100000.0;
               double l_mid;

               std::string bifur = "none";
               std::vector<double> errs = {0, 0, 0};
               std::vector<double> errs_old = {1, 1, 1};

               do
               {
                   err_old = err_new;
                   l_mid = (l_right  + l_left) / 2;
                   std::vector<double> ls = {{ l_left, l_mid, l_right }};

                  #pragma omp parallel for
                   for (std::vector<double>::size_type i = 0;
                        i < ls.size();
                        ++i)
                   {
                      dlib::matrix<double> m = par;
                      m(0, 0) = ls[i];
                      if (iter == 0)
                      {
                        errs[i] = loss(m);
                      }
                      else
                      {
                        errs_old = errs;
                        if (i != 1)
                        {
                          if ((i == 0 && bifur == "left" ) ||
                              (i == 2 && bifur == "right")) {
                            errs[i] = errs_old[1];
                          }
                        }
                        else {
                          errs[i] = loss(m);
                        }
                      }
                   }
                   if (errs[1] <  errs[2])
                   {
                     l_right = l_mid;
                      bifur = "right";
                   }
                   else if (errs[1] <  errs[0])
                   {
                     l_left  = l_mid;
                     bifur = "left";
                   }
                   else if (errs[1] == errs[2])
                   {
                     l_right = l_mid;
                     bifur = "right";
                   }
                   else
                   {
                     #ifdef USE_RCPPARMADILLO
                     Rprintf("Error estimating parameter %d using bifurcation!", idx);
                     #else
                     std::cerr << "Error estimating parameter " << idx << " using bifurcation!" << std::endl;
                     #endif
                     // return smallest value so far
                     return l_mid;
                   }
                   err_new = errs[1];
               }
               while(std::abs(err_new - err_old) > epsilon && ++iter < niter);
               //
               return l_mid;
           }
    };
}
#endif //NETREG_OPTIM_HPP
