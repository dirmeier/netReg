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

        /**
         * Returns a graph_penalized_linear_model_data object
         *
         * @param x the design matrix in column first order
         * @param y the response matrix in column first order
         * @param gx the adjacency matrix for variables of x
         * @param gy the adjacency matrix for variables of y
         * @param xdim the dimensionality of matrix x, where the first index
         *  is the number of rows and the second index the number of columns
         * @param ydim the dimensionality of matrix y, where the first index
         *  is the number of rows and the second index the number of columns
         * @param lambda regularization parameter for the LASSO.
         *  Only considered if do_lambda is false.
         * @param psigx regularization parameter for gx. Only considered if
         *  do_psigx is false.
         * @param psigy regularization parameter for gy. Only considered if
         *  do_psigy is false.
         * @param do_psigx if true uses regularization of GX and psigx. Otherwise no
         *  regularization will be used.
         * @param do_psigy if true uses regularization of GY and psigy. Otherwise no
         *  regularization will be used.
         * @param niter number of maximal iterations for coordinate descent
         * @param thresh threshold for convergence of likelihood
         * @param family the family of the likelihood
         *
         * @return returns a graph_penalized_linear_model_cv_data object
         */
        static graph_penalized_linear_model_data build_data(
          double* x, double* y,
          double* gx, double* gy,
          int* xdim, int* ydim,
          double lambda, double psigx, double psigy,
          bool do_psigx, bool do_psigy,
          int niter, double thresh,
          std::string& fam)
        {
            netreg::graph_penalized_linear_model_data data(
              x, y,
              gx, gy,
              xdim[0], xdim[1], ydim[1],
              lambda, psigx, psigy,
              do_psigx, do_psigy,
              niter, thresh,
              family(fam)
            );

            return data;
        }

        /**
         * Returns a graph_penalized_linear_model_cv_data object
         *
         * @param x the design matrix in column first order
         * @param y the response matrix in column first order
         * @param gx the adjacency matrix for variables of x
         * @param gy the adjacency matrix for variables of y
         * @param xdim the dimensionality of matrix x, where the first index
         *  is the number of rows and the second index the number of columns
         * @param ydim the dimensionality of matrix y, where the first index
         *  is the number of rows and the second index the number of columns
         * @param lambda regularization parameter for the LASSO.
         *  Only considered if do_lambda is false.
         * @param psigx regularization parameter for gx. Only considered if
         *  do_psigx is false.
         * @param psigy regularization parameter for gy. Only considered if
         *  do_psigy is false.
         * @param do_lambda boolean if the optimal regularization parameter
         *  for lambda should be estimated. If true disregards lambda and
         *  instead estimated.
         * @param do_psigx boolean if the optimal regularization parameter
         *  psigx should be estimated. If true disregards psigx and
         *  instead estimated.
         * @param do_psigy boolean if the optimal regularization parameter
         *  psigy should be estimated. If true disregards psigy and
         *  instead estimated.
         * @param niter number of maximal iterations for coordinate descent
         * @param thresh threshold for convergence of likelihood
         * @param lenfoldid length of the ptr foldids
         * @param foldids indexes of cross-validation folds
         * @param family the family of the likelihood
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
                  lambda, psigx, psigy,
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
                  lambda, psigx, psigy,
                  niter, thresh,
                  nfolds,
                  family(fam));

                return data;
            }
        }

    private:
        static const std::string GAUSSIAN;
        static const std::string BINOMIAL;
        static const std::string FAMILY_ERROR;

        static family family(std::string& fam)
        {
            enum family f =
              fam == GAUSSIAN ? family::GAUSSIAN :
              fam == BINOMIAL ? family::BINOMIAL :
              family::NONE;

            if (f == family::NONE)
            {
                #ifdef USE_RCPPARMADILLO
                Rcpp::stop(FAMILY_ERROR + "\n");
                #else
                throw std::invalid_argument(FAMILY_ERROR + "\n");
                #endif

            }

            return f;
        }

    };
};

#endif //NETREG_DATA_FACTORY_HPP

const std::string netreg::data_factory::GAUSSIAN = "gaussian";
const std::string netreg::data_factory::BINOMIAL = "binomial";
const std::string netreg::data_factory::FAMILY_ERROR =
  "Wrong family given. Choose one of " +
  netreg::data_factory::GAUSSIAN + "/" +
  netreg::data_factory::BINOMIAL;
