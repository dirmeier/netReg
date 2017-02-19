#include <iostream>

#include "family.hpp"
#include "edgenet_gaussian.hpp"
#include "graph_penalized_linear_model_data.hpp"

#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>

#include "armadillo"

int main()
{
    boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
    boost::normal_distribution<> nd(0.0, 1.0);
    boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<> > var_nor(rng, nd);

    const int n = 1000;
    const int p = 10000;
    const int q = 10;

    double *x = new double[n*p];
    double *y = new double[n*q];
    double *gx = new double[p*p];
    double *gy = new double[q*q];

    for (int j = 0; j < n*p; ++j) x[j] = var_nor();
    for (int j = 0; j < n*q; ++j) y[j] = var_nor();
    for (int k = 0; k < p*p; ++k) gx[k] = 1;
    for (int k = 0; k < q*q; ++k) gy[k] = 1;

    netreg::graph_penalized_linear_model_data data
        (x,y, gx, gy,
         n,p, q,
         10, 1.0,
         10, 10,
         100000, 0.00001, netreg::family::GAUSSIAN);

    netreg::edgenet_gaussian e;
    arma::Mat<double> m = e.run(data);

    std::cout << m << std::endl;


    delete [] x;
    delete [] y;
    delete [] gx;
    delete [] gy;

    return 0;

}
