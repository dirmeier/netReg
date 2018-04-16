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
#include "graph_functions.hpp"
#include "graph_model_data.hpp"
#include "graph_model_cv_data.hpp"
#include "family.hpp"


namespace netreg
{
    class data_factory
    {
    public:
        static graph_model_data build_data(
          double* x, double* y, double* gx, double* gy,
          int n, int p, int q, std::string& fam)
        {
            arma::Mat<double> xm(x, n, p, false, true);
            arma::Mat<double> ym(y, n, q, false, true);
            arma::Mat<double> gxm(gx, p, p, false, true);
            arma::Mat<double> gym(gy, q, q, false, true);

            netreg::graph_model_data data(
              xm, ym, gxm, gym,
              get_family(fam));

            return data;
        }

        static graph_model_data build_data(
          double* x, double* y, double* gx, double* gy,
          int* xdim, int* ydim, std::string& fam)
        {
            const int n = xdim[0];
            const int p = xdim[1];
            const int q = ydim[1];

            arma::Mat<double> xm(x, n, p, false, true);
            arma::Mat<double> ym(y, n, q, false, true);
            arma::Mat<double> gxm(gx, p, p, false, true);
            arma::Mat<double> gym(gy, q, q, false, true);

            netreg::graph_model_data data(
              xm, ym, gxm, gym,
              get_family(fam));

            return data;
        }

        static graph_model_cv_data build_cv_data(
          double* x, double* y, double* gx, double* gy,
          int n, int p, int q, std::string& fam, int nfolds)
        {
            arma::Mat<double> xm(x, n, p, false, true);
            arma::Mat<double> ym(y, n, q, false, true);
            arma::Mat<double> gxm(gx, p, p, false, true);
            arma::Mat<double> gym(gy, q, q, false, true);

            graph_model_cv_data data(
              xm, ym, gxm, gym, nfolds,
              get_family(fam));

            return data;
        }

        static graph_model_cv_data build_cv_data(
          double* x, double* y, double* gx, double* gy,
          int* xdim, int* ydim, std::string& fam,
          int nfolds, int lenfoldid, int* foldids)
        {
            const int n = xdim[0];
            const int p = xdim[1];
            const int q = ydim[1];

            arma::Mat<double> xm(x, n, p, false, true);
            arma::Mat<double> ym(y, n, q, false, true);
            arma::Mat<double> gxm(gx, p, p, false, true);
            arma::Mat<double> gym(gy, q, q, false, true);

            if (lenfoldid == xdim[0])
            {
                graph_model_cv_data data(
                  xm, ym, gxm, gym,
                  nfolds, foldids, get_family(fam));

                return data;
            }
            else
            {
                graph_model_cv_data data(
                  xm, ym, gxm, gym, nfolds, get_family(fam));

                return data;
            }
        }

    private:
        static const std::string GAUSSIAN;
        static const std::string BINOMIAL;
        static const std::string FAMILY_ERROR;

        static enum family get_family(std::string& fam)
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
