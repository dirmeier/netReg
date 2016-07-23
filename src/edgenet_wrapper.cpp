/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "edgenet_wrapper.hpp"

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "cv_set.hpp"
#include "graph_penalized_linear_model_data.hpp"
#include "edgenet.hpp"
#include "edgenet_model_selection.hpp"

void do_gauss_edgenet_(double *const X, double *const Y,
                       double *const GX, double *const GY,
                       const int n, const int p, const int q,
                       const double lambda,
                       const double psigx, const double psigy,
                       const int n_iter, const double thresh)
{
    netreg::graph_penalized_linear_model_data data(X, Y,
                                                   GX, GY,
                                                   n, p, q,
                                                   lambda, 1.0,
                                                   psigx, psigy,
                                                   n_iter, thresh);
    netreg::edgenet e;
    e.run(data);
    B_ = data.coefficients().begin();
    mu_ = data.intercept().begin();
}

void do_gauss_cv_edgenet_(double *const X, double *const Y,
                          double *const GX, double *const GY,
                          const int n, const int p, const int q,
                          const double psigx, const double psigy,
                          const int n_iter, const double thresh,
                          const int n_folds, int *const foldid,
                          const int n_foldid)
{
    netreg::graph_penalized_linear_model_data data(X, Y,
                                                   GX, GY,
                                                   n, p, q,
                                                   -1, 1.0,
                                                   psigx, psigy,
                                                   n_iter, thresh);
    netreg::edgenet_model_selection e;
    std::vector<double> pop;
    if (n_foldid == n)
        pop = e.regularization_path(data, foldid);
    else
        pop = e.regularization_path(data, n_folds);
    lamb_ = pop[0];
    psi_gx_ = pop[1];
    psi_gy_ = pop[2];
}
