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

#include "stat_functions.hpp"

namespace netreg
{
    cvector<double> intercept(matrix<double> &X, matrix<double> &Y,
                              matrix<double> &B)
    {
        matrix<double> terr = (Y - (X * B)).t();
        cvector<double> rep(Y.n_rows, arma::fill::ones);
        cvector<double> intr = (terr * rep) / Y.n_rows;
        return intr;
    }

    double partial_residual(matrix<double> &X,
                            matrix<double> &Y,
                            matrix<double> &B,
                            const int row,
                            const int pi,
                            const int qi)
    {
        double res = -X(row, pi) * B(pi, qi);
        const int length = static_cast<int>(B.n_rows);
        for (int p = 0; p < length; ++p)
            res += X(row, p) * B(p, qi);
        return Y(row, qi) - res;
    }

    double l_pls(matrix<double> &TXX,
                 matrix<double> &TXY,
                 matrix<double> &cfs,
                 const int cidx,
                 const int qi,
                 const int P)
    {

        double s = TXY(cidx, qi) + (TXX(cidx, cidx) * cfs(cidx, cidx));
        for (int j = 0; j < P; ++j)
        {
            // we only initialized the lower triangular matrix<double> of TXX
            if (j > cidx)
                s -= TXX(j, cidx) * cfs(j, cidx);
            else
                s -= TXX(cidx, j) * cfs(j, cidx);
        }
        return s;
    }

    double pls(matrix<double> &TXX,
               matrix<double> &TXY,
               matrix<double> &cfs,
               const int pi,
               const int qi,
               const int P,
               const bool lower)
    {

        if (lower)
            return l_pls(TXX, TXY, cfs, pi, qi, P);
        else
        {
            double s = TXY(pi, qi) + (TXX(pi, pi) * cfs(pi, qi))
                       - arma::accu(TXX.row(pi) * cfs.col(qi));
            return s;
        }
    }
}