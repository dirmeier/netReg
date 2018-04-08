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
#include "edgenet_binomial.hpp"

#include "math_functions.hpp"
#include "stat_functions.hpp"

namespace netreg
{
    SEXP edgenet_binomial::run(graph_penalized_linear_model_data &data) const
    {
        //        const int P = data.covariable_count();
        //        const int Q = data.response_count();
        //        matrix<double> &coef = data.coefficients();
        //        matrix<double> old_coef(P, Q);
        //        cvector<double> &intr = data.intercept();
        //        const double thresh = data.threshold();
        //        const int niter = data.max_iter();
        //        const double lambda = data.lambda();
        //        const double alpha = data.alpha();
        //        const double psigx=data.psigx();
        //        const double psigy=data.psigy();
        //        matrix<double> &TXX= data.txx();
        //        matrix<double> &TXY = data.txy();
        //        matrix<double> &LX = data.lx();
        //        matrix<double> &LY = data.ly();
        //        // calculate intercepts of the linear model
        //        intr = intercept(data.design(), data.response(), coef);
        return R_NilValue;
    }

    arma::Mat<double> edgenet_binomial::run_cv(
      graph_penalized_linear_model_cv_data &data,
      const double lambda,
      const double psigx,
      const double psigy,
      cv_fold &fold) const
    {
        arma::Mat<double> coef(10, 10, arma::fill::ones);
        // TODO
        return coef;
    }
}
