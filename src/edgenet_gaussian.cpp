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
    arma::Mat<double> edgenet_gaussian::run(
      graph_penalized_linear_model_data& data) const
    {
        // load shrinkage coefficients
        const double lambda = data.lambda();
        const double alpha  = data.alpha();
        const double psigx  = data.psigx();
        const double psigy  = data.psigy();

        // load square matrices
        std::vector<arma::rowvec>& txx_rows = data.txx_rows();
        arma::Mat<double>& txy              = data.txy();

        return mccd_(data, lambda, alpha, psigx, psigy, txx_rows, txy);
    }

    arma::Mat<double> edgenet_gaussian::run_cv(
      graph_penalized_linear_model_cv_data& data,
      const double lambda,
      const double alpha,
      const double psigx,
      const double psigy,
      cv_fold& fold) const
    {
        // load square matrices
        // but load the TRAINING matrices of the fold
        std::vector<arma::rowvec>& train_txx_rows = fold.train_txx_rows();
        arma::Mat<double>& train_txy              = fold.train_txy();

        return mccd_(
          data, lambda, alpha, psigx, psigy, train_txx_rows, train_txy);
    }

    arma::Mat<double> edgenet_gaussian::mccd_(
      graph_penalized_linear_model_data& data,
      const double lambda,
      const double alpha,
      const double psigx,
      const double psigy,
      std::vector<arma::rowvec>& txx_rows,
      arma::Mat<double>& txy) const
    {
        /*
         * Load variables from data-set
         */
        const int P         = data.covariable_count();
        const int Q         = data.response_count();
        const double thresh = data.threshold();
        const int niter     = data.max_iter();

        // load graph laplacians
        std::vector<arma::Row<double>>& lx_rows = data.lx_rows();
        arma::Mat<double>& ly                   = data.ly();

        // setup coefficient matrix
        arma::Mat<double> coef(P, Q, arma::fill::ones);
        arma::Mat<double> old_coef(P, Q);
        std::vector<arma::rowvec> coef_rows(static_cast<unsigned int>(P));

        // TODO: method for this
        // safe an extra set of rowvectors so that access is faster
        #pragma omp parallel for
        for (std::vector<arma::Row<double>>::size_type i = 0;
             i < coef.n_rows;
             ++i)
        {
            coef_rows[i] = coef.row(i);
        }

        // TODO methods for this
        // weighted penalization param of Elastic-net
        const double lalph = alpha * lambda;
        // normalization for soft-thresholding
        const double enorm = 1.0 + lambda * (1 - alpha);

        for (int qi = 0; qi < Q; ++qi)
        {
            uccd_(P, Q,
                  qi,
                  thresh, niter,
                  lalph, enorm,
                  psigx, psigy,
                  txx_rows, txy,
                  lx_rows, ly,
                  coef, old_coef, coef_rows);
#ifdef USE_RCPPARMADILLO
            if (qi % 100 == 0)
            {
                Rcpp::checkUserInterrupt();
            }
#endif
        }

        return coef;
    }

    void edgenet_gaussian::uccd_(const int P,
                                 const int Q,
                                 const int qi,
                                 const double thresh,
                                 const int niter,
                                 const double lalph,
                                 const double enorm,
                                 const double psigx,
                                 const double psigy,
                                 std::vector<arma::rowvec>& txx_rows,
                                 arma::Mat<double>& txy,
                                 std::vector<arma::rowvec>& lx_rows,
                                 arma::Mat<double>& ly,
                                 arma::Mat<double>& coef,
                                 arma::Mat<double>& old_coef,
                                 std::vector<arma::rowvec>& coef_rows) const
    {
        // iteration counter
        int iter    = 0;
        double s    = 0.0;
        double norm = 0.0;
        // do while estimation of params does not converge
        do
        {
            // fix bnew_i and calculate least-squares
            // coefficient on partial residual
            for (int pi = 0; pi < P; ++pi)
            {
                // TODO methodize
                // safe current estimate of coefficients
                old_coef(pi, qi) = coef(pi, qi);
                // TODO: no void stuff :(
                set_params(s, norm,
                           P, Q,
                           pi, qi,
                           psigx, psigy,
                           txx_rows[pi], txy,
                           lx_rows[pi], ly,
                           coef, coef_rows[pi]);
                // soft-thresholded version of estimate
                const double d = softnorm(s, lalph, enorm * norm);
                // TODO: METHOD for this
                coef(pi, qi) = d;
                coef_rows[pi](qi) = d;
#ifdef USE_RCPPARMADILLO
                if (iter % 100 == 0)
                {
                    Rcpp::checkUserInterrupt();
                }
#endif
            }
        }
        // TODO method for accumulation
        while (arma::accu(arma::abs(coef.col(qi) - old_coef.col(qi))) >
                   thresh && iter++ < niter);
    }
}
