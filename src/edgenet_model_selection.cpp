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

#include "edgenet_gaussian_model_selection.hpp"

#include "family.hpp"
#include "optim.hpp"
#include "edgenet_gaussian_loss_function.hpp"

namespace netreg
{
    std::map<std::string, double> model_selection(
       graph_model_cv_data& cv_data, params& pars)
    {
        optim opt;

        std::vector<double> start{0, 0, 0};
        std::vector<double> lower_bound{0.0, 0.0, 0.0};
        std::vector<double> upper_bound{100.0, 10000.0, 10000.0};
        const double rad_start = 0.49, rad_end = pars.optim_epsilon();

        std::map<std::string, double> res;
        switch (cv_data.data().distribution_family())
        {
            case family::GAUSSIAN:
            default:
            {
                return opt.bobyqa<edgenet_gaussian_loss_function>(
                  cv_data,
                  pars,
                  start,
                  lower_bound,
                  upper_bound,
                  rad_start,
                  rad_end,
                  pars.optim_niter());
            }
        }

        return res;
    }
}
