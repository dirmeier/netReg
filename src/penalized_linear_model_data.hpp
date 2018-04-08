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


#ifndef NETREG_PENALIZED_LINEAR_MODEL_DATA_HPP
#define NETREG_PENALIZED_LINEAR_MODEL_DATA_HPP

#include "linear_model_data.hpp"

#include <string>
#include "family.hpp"

namespace netreg
{
    /**
     * Data-structure for all required data for penalized regression.
     */
    class penalized_linear_model_data: public linear_model_data
    {
    public:
        penalized_linear_model_data()
        {}

        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param n the number of samples (nrows X/Y)
         * @param p the number of covariables (ncols X)
         * @param q the number of responses (ncol Y)
         * @param niter max number of iterations in case estimation of the
         * coefficients does not converge
         * @param thresh convergence threshold
         * @param is_column_first boolean if given matrices are column-first or
         * not
         */
        penalized_linear_model_data(
          double *const x, double *const y, int n, int p, int q,
          double lambda, int niter, double thresh,
          const enum family fam):
          linear_model_data(x, y, n, p, q, niter, thresh, fam),
          LAMBDA(lambda)
        {
        }

        /**
         * Getter for the penalization parameter of the lasso.
         *
         * @return the penalization term for the LASSO
         */
        const double lambda()
        {
            return LAMBDA;
        }

    protected:
        double LAMBDA;  // penalization term for LASSO
    };
}
#endif  // NETREG_PENALIZEDLINEARMODELDATA_HPP
