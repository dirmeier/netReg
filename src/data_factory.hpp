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
         * Creates a graph_model_data object
         *
         * @param x the design matrix in column first order
         * @param y the response matrix in column first order
         * @param gx the adjacency matrix for variables of x
         * @param gy the adjacency matrix for variables of y
         * @param xdim the dimensionality of matrix x, where the first index
         *  is the number of rows and the second index the number of columns
         * @param ydim the dimensionality of matrix y, where the first index
         *  is the number of rows and the second index the number of columns
         * @param family the family of the likelihood
         *
         * @return returns a graph_model_data object
         */
        static graph_model_data build_data(
          double* x, double* y, double* gx, double* gy,
          int* xdim, int* ydim, std::string& fam)
        {
            const int n = xdim[0];
            const int p = xdim[1];
            const int q = ydim[1];

            netreg::graph_model_data data(
              arma::Mat<double>(x, n, p, false, true),
              arma::Mat<double>(y, n, q, false, true),
              arma::Mat<double>(gx, p, p, false, true),
              arma::Mat<double>(gy, q, q, false, true),
              family(fam));

            return data;
        }

        /**
         * Creats a graph_model_cv_data object
         *
         * @param x the design matrix in column first order
         * @param y the response matrix in column first order
         * @param gx the adjacency matrix for variables of x
         * @param gy the adjacency matrix for variables of y
         * @param xdim the dimensionality of matrix x, where the first index
         *  is the number of rows and the second index the number of columns
         * @param ydim the dimensionality of matrix y, where the first index
         *  is the number of rows and the second index the number of columns
         * @param lenfoldid length of the ptr foldids
         * @param foldids indexes of cross-validation folds
         * @param family the family of the likelihood
         *
         * @return returns a graph_model_cv_data object
         */
        static graph_model_cv_data build_cv_data(
          double* x, double* y, double* gx, double* gy,
          int* xdim, int* ydim, std::string& fam,
          int nfolds, int lenfoldid, int* foldids)
        {

            const int n = xdim[0];
            const int p = xdim[1];
            const int q = ydim[1];

            if (lenfoldid == xdim[0])
            {
                graph_model_cv_data data(
                  arma::Mat<double>(x, n, p, false, true),
                  arma::Mat<double>(y, n, q, false, true),
                  arma::Mat<double>(gx, p, p, false, true),
                  arma::Mat<double>(gy, q, q, false, true),
                  foldids,
                  family(fam));

                return data;
            }
            else
            {
                graph_model_cv_data data(
                  arma::Mat<double>(x, n, p, false, true),
                  arma::Mat<double>(y, n, q, false, true),
                  arma::Mat<double>(gx, p, p, false, true),
                  arma::Mat<double>(gy, q, q, false, true),
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
