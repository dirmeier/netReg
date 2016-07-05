/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "edgenet_wrapper.hpp"

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "graph_penalized_linear_model_data.hpp"
#include "edgenet.hpp"
#include "edgenet_model_selection.hpp"

void do_gauss_edgenet_(double *const X, double *const Y,
         double *const GX, double *const GY,
         const int N, const int P, const int Q,
         const double LAMBDA, const double PSI_GX, const double PSI_GY,
         const int N_ITER, const double THRESH)
{
    netreg::graph_penalized_linear_model_data data(X, Y,
                                                   GX, GY,
                                                   N, P, Q,
                                                   LAMBDA, 1.0, PSI_GX,
                                                   PSI_GY,
                                                   N_ITER, THRESH);
    // create an edgenet object and fit the model
    netreg::edgenet e;
    e.run(data);
    B_  = data.coefficients().begin();
    mu_ = data.intercept().begin();
}

void do_gauss_cv_edgenet_(double *const X, double *const Y,
                          double *const GX, double *const GY,
                          const int N, const int P, const int Q,
                          const double LAMBDA, const double PSI_GX, const double PSI_GY,
                          const int N_ITER, const double THRESH,
                          const int N_FOLDS, int * const foldid)
{

}
