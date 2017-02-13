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


#include "edgenet.hpp"
#include "edgenet_gaussian.hpp"

namespace netreg
{
    SEXP edgenet::run(graph_penalized_linear_model_data &data) const
    {
        netreg::edgenet_gaussian edge;
        return edge.run(data);
    }

    arma::Mat<double> edgenet::run_cv
        (graph_penalized_linear_model_cv_data &data,
         const double lambda, const double alpha,
         const double psigx,  const double psigy,
         cv_fold &fold) const
    {
        netreg::edgenet_gaussian edge;
        return edge.run_cv(data, lambda, alpha, psigx, psigy, fold);
    }
}
