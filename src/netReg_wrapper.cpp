/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include <Rcpp.h>
#include "graph_penalized_linear_model_data.hpp"

/**
 * Implementation of Edgenet, a edge-based regularized regression model.
 *
 * @param XS the (ns x ps)-dimensional design matrix
 * @param YS the (ns x qs)-dimensional response matrix
 * @param GXS the (ps x ps)-prior graph for XS
 * @param GYS the (qs x qs)-prior graph for YS
 * @param ns number of observations/samples/rows
 * @param ps number of covariates
 * @param qs number of responses
 * @param lamdbass penalization value for LASSO
 * @param psi_gxs weighting value of GX
 * @param psi_gys weighting value of GY
 * @param niters max number of iterations if parameter estimation does not converge in time
 * @param threshs convergence threshold
 */
// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::export(name=".edgenet.cpp")]]
Rcpp::List edgenet_rcpp_(
    const Rcpp::NumericMatrix& XS, const Rcpp::NumericMatrix& YS,
    const Rcpp::NumericMatrix& GXS, const Rcpp::NumericMatrix& GYS,
    const int n, const int p, const int q,
    const double lambda, const double psigx, const double psigy,
    const int n_iter, const double thresh,
    const Rcpp::CharacterVector& familys)
{
    double * X = REAL(XS);
    double * Y = REAL(YS);
    double * GX = REAL(GXS);
    double * GY = REAL(GYS);
    netreg::graph_penalized_linear_model_data data(X, Y,
                                                   GX, GY,
                                                   n, p, q,
                                                   lambda, 1.0,
                                                   psigx, psigy,
                                                   n_iter, thresh);
    return Rcpp::List::create(1);
}

