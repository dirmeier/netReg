/**
 * Author: Simon Dirmeier
 * Date: 8/18/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#include <iostream>
int main()
{
    arma::Mat<double> a(10, 10, arma::fill::zeros);
    arma::Mat<double> b = 10 - a;
    for (unsigned int i = 0; i < b.n_rows; ++i)
    {
        for (unsigned int j = 0; j < b.n_cols; ++j)
            std::cout << b(i,j) << " ";
        std::cout << "\n";
    }
}