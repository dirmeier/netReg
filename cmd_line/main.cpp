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

    const int n = 100;
    const int p = 100;
    const int q = 10;

    double *x = new double[n * p];
    double *y = new double[n * q];
    double *gx = new double[p * p];
    double *gy = new double[q * q];

    for (int j = 0; j < n * p; ++j) x[j] = 1;
    for (int j = 0; j < n * q; ++j) y[j] = 2;
    for (int k = 0; k < p * p; ++k) gx[k] = 1;
    for (int k = 0; k < q * q; ++k) gy[k] = 1;

  /*  netreg::graph_penalized_linear_model_data dat
        (x, y, gx, gy, n, p, q,
         10, 1.0, 1.0, 1.0, 100000, 0.0001,
         netreg::family::GAUSSIAN);

    netreg::edgenet_gaussian ed;
    ed.run(dat);
    std::cout << "done1" << std::endl;*/

    netreg::graph_penalized_linear_model_cv_data data
        (x, y, gx, gy, n, p, q,
         -1, 1.0, 0.0, 0.0, 100000, 0.0000000001, 5,
         netreg::family::GAUSSIAN);

    netreg::edgenet_gaussian_model_selection e;

    std::map<std::string, double> m = e.regularization_path(data, 1000, 0.001);
    std::cout << "Approx: " << m["lambda"] << " " << m["psigx"] << " " << m["psigy"] << std::endl;

    // std::map<std::string, double> m = e.regularization_path(data, false, 1000, 0.001);
    // std::cout << "True: " <<  m["lambda"] << " " << m["psigx"] << " " << m["psigy"] << std::endl;

    delete[] x;
    delete[] y;
    delete[] gx;
    delete[] gy;


    return 0;

}
