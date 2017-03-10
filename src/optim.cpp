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

#include <cmath>;

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#include <iostream>
#endif

#include "../inst/include/dlib/optimization.h"

#include "graph_penalized_linear_model_cv_data.hpp"
#include "cv_set.hpp"
#include "edgenet_gaussian_loss_function.hpp"

namespace netreg
{

        template<typename loss_function>
        std::map<std::string, double> optim::bobyqa
            (graph_penalized_linear_model_cv_data &data,
             std::vector<double> &start,
             std::vector<double> &lower_bound,
             std::vector<double> &upper_bound,
             const double radius_start,
             const double radius_stop,
             const int niter)
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

        template<typename loss_function>
        std::map<std::string, double> optim::bifurcation
          (graph_penalized_linear_model_cv_data &data,
           std::vector<double> &lower_bound,
           std::vector<double> &upper_bound,
           const double epsilon,
           const int niter)
           {
             const int sz = static_cast<int>(start.size());
             loss_function loss(data);
             // convert to dlib objects
             dlib::matrix<double> par(sz, 1);
             for (int i = 0; i < sz; ++i) par(i, 0) = start[i];
             // minimize the loss_function
             try
             {
                for (std::vector<double>::size_type i = 0;
                     i < upper_bound.size();
                     ++i)
                  {
                      double opt = blockwise_bifurcation<loss_function>(
                        par, i, liss,
                        lower_bound, upper_bound,
                        epsilon, niter);
                      par(i, 0) = opt;
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

           template<typename loss_function>
           double optim::blockwise_bifurcation
             (const dlib::matrix<double>& par, const int idx,
              loss_function &loss,
              std::vector<double> &lower_bound,
              std::vector<double> &upper_bound,
              const double epsilon,
              const int niter)
              {

                  int iter = 0;
                  double l_left  = lower_bound(i, 1);
                  double l_right = upper_bound(i, 1);
                  double err_old = 100000.0;
                  double err_new = 0;

                  do
                  {
                      err_old = err_new;

                      double l_mid = (l_right  - l_left) / 2;

                      std::vector<double> ls = {{ l_left, l_mid, l_right }};
                      std::vector<double> errs(3);

                      // parallelize
                      for (std::vector<double>::size_type i = 0;
                           i < ls.size();
                           ++i)
                      {
                          dlib::matrix<double> m = par;
                          m(i, 0) = ls[i];
                          errs[i] = loss(m);
                      }

                      if      (errs[1] <  errs[2])  l_right = l_mid;
                      else if (errs[1] <  errs[0])  l_left  = l_mid;
                      else if (errs[1] == errs[2]]) l_right = l_mid;
                      else
                      {
                        #ifdef USE_RCPPARMADILLO
                        Rprintf("Error estimating parameter %d using bifurcation!", i);
                        #else
                        std::cerr << "Error estimating parameter " << i << " using bifurcation!" << std::endl;
                        #endif
                        // return smallest value so far
                        return l_mid;
                      }

                      err_new = errs[1];
                  }
                  while(std::abs(err - err_old) > epsilon && ++iter < niter);

                  return l_mid;
              }
    };
}
