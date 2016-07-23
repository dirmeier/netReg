/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "stat_functions.hpp"

#include <iostream>

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
//        std::cout << "lower" << lower << std::endl;
//        std::cout << "pi" << pi << std::endl;
//        std::cout << "qi" << qi << std::endl;
//        std::cout << "P" << P << std::endl;
//        for (int i = 0; i < TXX.n_rows; ++i)
//        {
//            for (int j = 0; j < TXX.n_cols; ++j)
//            {
//                std::cout << TXX(i, j) << " ";
//            }
//            std::cout <<  "\n";
//        }
//        std::cout <<  "\n";
//        for (int i = 0; i < TXY.n_rows; ++i)
//        {
//            for (int j = 0; j < TXY.n_cols; ++j)
//            {
//                std::cout << TXY(i, j) << " ";
//            }
//            std::cout <<  "\n";
//        }
//        std::cout <<  "\n";

        if (lower)
            return l_pls(TXX, TXY, cfs, pi, qi, P);
        else
        {
            double s = TXY(pi, qi) + (TXX(pi, pi) * cfs(pi, qi));
            for (int j = 0; j < P; j++)
                s -= (TXX)(pi, j) * cfs(j, qi);
            return s;
        }
    }
}