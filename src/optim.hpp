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

#include "graph_penalized_linear_model_cv_data.hpp"

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
             const int niter);

        // todo
        template<typename loss_function>
        std::map<std::string, double> bifurcation
          (graph_penalized_linear_model_cv_data &data,
           std::vector<double> &start,
           std::vector<double> &lower_bound,
           std::vector<double> &upper_bound,
           const double epsilon,
           const int niter);
    };
}
#endif //NETREG_OPTIM_HPP
