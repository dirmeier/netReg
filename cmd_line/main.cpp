#include <iostream>
#include <map>
#include <string>

#include "../src/family.hpp"
#include "../src/edgenet_gaussian.hpp"
#include "../src/graph_penalized_linear_model_cv_data.hpp"
#include "../src/edgenet_gaussian_model_selection.hpp"

#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>


int main()
{
    boost::mt19937 rng; // I don't seed it on purpouse (it's not relevant)
    boost::normal_distribution<> nd(0.0, 1.0);
    boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<> > var_nor(rng, nd);

    const int n = 10;
    const int p = 100;
    const int q = 1;

    double *x = new double[n*p];
    double *y = new double[n*q];
    double *gx = new double[p*p];
    double *gy = new double[q*q];

    for (int j = 0; j < n*p; ++j) x[j] = var_nor();
    for (int j = 0; j < n*q; ++j) y[j] = var_nor();
    for (int k = 0; k < p*p; ++k) gx[k] = 1;
    for (int k = 0; k < q*q; ++k) gy[k] = 1;

    netreg::graph_penalized_linear_model_cv_data data
        (x, y, gx, gy, n, p, q,
         -1, 1.0, -1, 0, 100000, 0.0000000001, 5,
         netreg::family::GAUSSIAN);

    netreg::edgenet_gaussian_model_selection e;
    std::map<std::string, double> m = e.regularization_path(data);

    delete [] x;
    delete [] y;
    delete [] gx;
    delete [] gy;

    return 0;

}
