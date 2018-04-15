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

#ifndef NETREG_MODEL_DATA_HPP
#define NETREG_MODEL_DATA_HPP

#include <vector>
#include <string>
#include <utility>

#ifdef USE_RCPPARMADILLO
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#else
#include "armadillo"
#endif
#include "family.hpp"

namespace netreg
{
    /**
     * Data-structure for all required data for a model
     */
    class model_data
    {
    public:

        /**
         * Protected constructor in order to avoid instantiation.
         *
         * @param x the design matrix
         * @param y the response matrix
         * @param fam of distribution of y
         */
        model_data(
          arma::Mat<double>& x, arma::Mat<double>& y, const family fam):
          N(x.n_rows), P(x.n_cols), Q(y.n_cols),
          X(x.memptr(), N, P, false, true),
          Y(y.memptr(), N, Q, false, true),
          TXY(P, Q),
          txx_rows_(P),
          family_(fam)
        {
            arma::Mat<double> TX = X.t();
            arma::Mat<double> TXX = TX * X;
            TXY = TX * Y;

            for (std::vector<arma::Row<double> >::size_type i = 0;
                 i < TXX.n_rows;
                 ++i)
            {
                txx_rows_[i] = TXX.row(i);
            }
        }

        /**
         * Getter for the family of the distribution of the response.
         *
         * @return returns the family
         */
        family distribution_family()
        {
            return family_;
        }

        /**
         * Getter for number of samples.
         *
         * @return returns the number of samples
         */
        int sample_count()
        {
            return N;
        }

        /**
         * Getter for number of responses.
         *
         * @return returns the number of responses
         */
        int response_count()
        {
            return Q;
        }

        /**
         * Getter for number of covariables.
         *
         * @return returns the number of covariables
         */
        int covariable_count()
        {
            return P;
        }

        arma::rowvec& txx_row(const int i)
        {
            return txx_rows_[i];
        }

        std::vector<arma::rowvec>& txx_rows()
        {
            return txx_rows_;
        }

        /**
         * Getter for the design matrix.
         *
         * @return a reference to the design matrix
         */
        arma::Mat<double>& design()
        {
            return X;
        }

        /**
         * Getter for the response matrix.
         *
         * @return a reference to the response matrix
         */
        arma::Mat<double>& response()
        {
            return Y;
        }

        /**
         * Getter for X'Y matrix.
         *
         * @return a reference to X'Y matrix.
         */
        arma::Mat<double>& txy()
        {
            return TXY;
        }

    protected:
        int N;                  // number of samples: n
        int P;                  // number of covariables: p
        int Q;                  // number of responses: q
        arma::Mat<double> X;    // (n x p)-dimensional design matrix
        arma::Mat<double> Y;    // (n x q)-dimensional response matrix
        arma::Mat<double> TXY;  // (p x q)-dimensional matrix: X'Y
        std::vector<arma::Row<double> > txx_rows_; // TXX as vector of rows
        enum family family_;    // family of distribution of y
    };
}
#endif
