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
#include "graph_functions.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <vector>
#include <cmath>

namespace netreg
{
    arma::Mat<double> laplacian(const double * x, const int n, const int m, const double px)
    {

        if (px == 0)
        {
            arma::Mat<double> lap(1, 1, arma::fill::zeros);
            return lap;
        }

        std::vector<double> degrees(n);
        #pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            double rowSum = 0.0;
            for (int j = 0; j < m; ++j)
                rowSum += x[i + n * j];
            degrees[i] = rowSum;
        }
        double *  laplacian = new double[n * m];
        // calculate normalized Laplacian matrix of source
        #pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                if (i == j && degrees[i] != 0)
                    laplacian[i + n * j] = 1 - (x[i + n * j] / degrees[i]);
                else if (i != j && x[i + n * j] != 0)
                    laplacian[i + n * j] =
                        -x[i + n * j] / std::sqrt(degrees[i] * degrees[j]);
                else
                    laplacian[i + n * j] = 0.0;
            }
        }
        arma::Mat<double> lap(laplacian, n, m, true, true);
        delete [] laplacian;
        return lap;
    }
}