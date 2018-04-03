/**
 * Author: Simon Dirmeier
 * Date: 02.04.18
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_DATA_FACTORY_HPP
#define NETREG_DATA_FACTORY_HPP

#ifdef USE_RCPPARMADILLO
#include <RcppArmadillo.h>
#endif

#include <string>
#include <stdexcept>
#include "graph_penalized_linear_model_data.hpp"
#include "graph_penalized_linear_model_cv_data.hpp"
#include "family.hpp"

namespace netreg
{
    class data_factory
    {
    public:

        static graph_penalized_linear_model_data build_model_data(
          double* x, double* y,
          double* gx, double* gy,
          int* xdim, int* ydim,
          double lambda, double psigx, double psigy,
          int niter, double thresh,
          std::string& fam)
        {
            netreg::graph_penalized_linear_model_data data(
              x, y,
              gx, gy,
              xdim[0], xdim[1], ydim[1],
              lambda, 1.0, psigx, psigy,
              niter, thresh,
              family(fam)
            );

            return data;
        }

        /**
         *  Returns a graph_penalized_linear_model_cv_data object
         *
         * @param x the design matrix in column first order
         * @param y the response matrix in column first order
         * @param gx
         * @param gy
         * @param xdim
         * @param ydim
         * @param lambda
         * @param psigx
         * @param psigy
         * @param do_lambda
         * @param do_psigx
         * @param do_psigy
         * @param niter
         * @param thresh
         * @param lenfoldid
         * @param foldids
         * @param family
         *
         * @return returns a graph_penalized_linear_model_cv_data object
         */
        static graph_penalized_linear_model_cv_data build_cv_data(
          double* x, double* y,
          double* gx, double* gy,
          int* xdim, int* ydim,
          double lambda, double psigx, double psigy,
          bool do_lambda, bool do_psigx, bool do_psigy,
          int niter, double thresh,
          int nfolds, int lenfoldid, int* foldids,
          std::string& fam)
        {
            if (lenfoldid == xdim[0])
            {
                graph_penalized_linear_model_cv_data data(
                  x, y, gx, gy,
                  xdim[0], xdim[1], ydim[1],
                  lambda, 1.0,
                  psigx, psigy,
                  niter, thresh,
                  foldids,
                  family(fam));

                return data;
            }
            else
            {
                graph_penalized_linear_model_cv_data data(
                  x, y, gx, gy,
                  xdim[0], xdim[1], ydim[1],
                  lambda, 1.0,
                  psigx, psigy,
                  niter, thresh,
                  nfolds,
                  family(fam));

                return data;
            }
        }

    private:
        family family(std::string& family)
        {
            family f =
              fam == "gaussian" ? family::GAUSSIAN :
              fam == "binomial" ? family::BINOMIAL :
              family::NONE;

            if (f == family::NONE)
            {
                #ifdef USE_RCPPARMADILLO
                Rcpp::stop("Wrong family given\n");
                #else
                throw std::invalid_argument(
                  "Error estimating optim shrinkage parameters.");
                #endif

            }

            return f;
        }

    };
};

#endif //NETREG_DATA_FACTORY_HPP
