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

#include "edgenet_model_selection.hpp"

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include <numeric>
#include <vector>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "family.hpp"
#include "edgenet_gaussian_loss_function.hpp"
#include "edgenet_binomial_loss_function.hpp"
#include "optim.hpp"

namespace netreg
{

    std::vector<double>
    edgenet_model_selection::regularization_path
        (graph_penalized_linear_model_cv_data &data)
    {
        optim opt;
        std::vector<double> start{0, 0, 0};
        std::vector<double> lower_bound{0.0, 0.0, 0.0};
        std::vector<double> upper_bound{100.0, 10000.0, 10000.0};
        const double rad_start = 0.49, rad_end = 1e-6;
        const int niter = 1000;
        switch (data.family())
        {
          case family::BINOMIAL:
                return opt.bobyqa<edgenet_binomial_loss_function>
                              (data, start, lower_bound, upper_bound,
                               rad_start, rad_end, niter);
          case family::GAUSSIAN:
          default:
                return opt.bobyqa<edgenet_gaussian_loss_function>
                              (data, start, lower_bound, upper_bound,
                               rad_start, rad_end, niter);

        }
    }
}