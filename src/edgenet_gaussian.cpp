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

#include "edgenet_gaussian.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "cv_set.hpp"

namespace netreg
{
    arma::Mat<double> edgenet_gaussian::run()
    {
        // load square matrices
        std::vector<arma::rowvec>& txx_rows = DATA_.txx_rows();
        arma::Mat<double>& txy = DATA_.txy();

        return mccd_(txx_rows, txy);
    }

    arma::Mat<double> edgenet_gaussian::run_cv(
      const double lambda,
      const double psigx,
      const double psigy,
      cv_fold& fold) const
    {
        // load the TRAINING squares matrices of the fold
        std::vector<arma::rowvec>& train_txx_rows = fold.train_txx_rows();
        arma::Mat<double>& train_txy = fold.train_txy();

        // set parameters we want to use for cross-validation
        set_lambda(lambda);
        set_psigx(psigx);
        set_psigy(psigy);

        return mccd_(train_txx_rows, train_txy);
    }

    arma::Mat<double> edgenet_gaussian::mccd_(
      std::vector<arma::rowvec>& txx_rows,
      arma::Mat<double>& txy) const
    {

        // setup coefficient matrix
        arma::Mat<double> B(P_, Q_, arma::fill::ones);
        arma::Mat<double> B_old(P_, Q_);
        std::vector<arma::rowvec> B_rows(static_cast<unsigned int>(P_));

        // safe an extra set of rowvectors so that access is faster
        for (std::vector<arma::Row<double>>::size_type i = 0;
             i < B.n_rows;
             ++i)
        {
            B_rows[i] = B.row(i);
        }

        for (unsigned int qi = 0; qi < Q_; ++qi)
        {
            uccd_(qi, txx_rows, txy, B, B_old, B_rows);

            #ifdef USE_RCPPARMADILLO
            if (qi % 100 == 0)
            {
                Rcpp::checkUserInterrupt();
            }
            #endif
        }

        return B;
    }

    void edgenet_gaussian::uccd_(const int qi,
                                 std::vector<arma::rowvec>& txx_rows,
                                 arma::Mat<double>& txy,
                                 arma::Mat<double>& B,
                                 arma::Mat<double>& B_old,
                                 std::vector<arma::rowvec>& B_rows) const
    {
        // iteration counter
        int iter = 0;
        double s = 0.0;
        double norm = 0.0;

        // do while estimation of params does not converge
        do
        {
            // fix bnew_i and calculate least-squares
            // coefficient on partial residual
            for (unsigned int pi = 0; pi < P_; ++pi)
            {
                // safe current estimate of coefficients
                B_old(pi, qi) = B(pi, qi);

                // TODO: no void stuff :(
                set_params(s, norm, pi, qi, txx_rows[pi], txy, B, B_rows[pi]);

                // soft-thresholded version of estimate
                const double d = softnorm(s, lambda_, norm);

                // TODO: METHOD for this
                B(pi, qi) = d;
                B_rows[pi](qi) = d;

                #ifdef USE_RCPPARMADILLO
                if (iter % 100 == 0)
                {
                    Rcpp::checkUserInterrupt();
                }
                #endif
            }
        }
        while (!converged_(B.col(qi), B_old.col(qi), iter++));
    }
}
