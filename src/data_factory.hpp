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
