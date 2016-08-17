/**
 * Author: Simon Dirmeier
 * Date: 06/08/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */


#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#include <iostream>
#include <cstdlib>

#include <boost/random.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include "../cv_set.hpp"
#include "../graph_penalized_linear_model_data.hpp"
#include "../gaussian_edgenet.hpp"
#include "../edgenet_model_selection.hpp"


// Version that fills a vector
template<class T>
void generate(T &generator, double * res, int sz)
{
    for(size_t i = 0; i < sz; ++i)
        res[i] = generator();
}

void generate( double * res,unsigned int sz)
{
    for(size_t i = 0; i < sz; ++i)
        res[i] = 1.0;
}

void norm(double * X, unsigned int n, unsigned int m)
{
    for(size_t i = 0; i < n; ++i)
        for(size_t j = 0; j < m; ++j)
            if (i == j)
                X[i * m + j] = 0.0;
            else
                X[i * m + j] = X[j * m + i];
}

int main(int argc, char** argv)
{

    boost::mt19937 rng(23);
    boost::normal_distribution<> snd(0.0, 1.0);
    boost::exponential_distribution<> ed(1.0);
    boost::variate_generator<boost::mt19937&,
        boost::normal_distribution<> > var_nor(rng, snd);
    boost::variate_generator<boost::mt19937&,
        boost::exponential_distribution<> > var_exp(rng, ed);

    const unsigned int n = atoi(argv[1]);
    const unsigned int p = atoi(argv[2]);
    const unsigned int q = atoi(argv[3]);
    double * X = new double[n * p];
    double * Y = new double[n * q];
    double * GX = new double[p * p];
    double * GY = new double[q * q];
    generate(var_nor, X, n*p);
    generate(var_nor, Y, n*q);
    generate(var_exp, GX, p*p);
    generate(var_exp, GY, q*q);
    norm(GX, p, p);
    norm(GY, q, q);
    netreg::graph_penalized_linear_model_data data(X, Y,
                                                   GX, GY,
                                                   n, p, q,
                                                   100, 1.0,
                                                   0.0, 0.0,
                                                   10000, 0.00001);
    netreg::edgenet e;
    e.run(data);
    delete [] X;
    delete [] GX;
    delete [] Y;
    delete [] GY;
    return 0;
}