/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include <Rcpp.h>
#include <string>

#include "graph_penalized_linear_model_data.hpp"
#include "edgenet.hpp"
#include "gaussian_edgenet.hpp"
#include "binomial_edgenet.hpp"
#include "edgenet_model_selection.hpp"

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
    const Rcpp::NumericMatrix& XS,
    const Rcpp::NumericMatrix& YS,
    const Rcpp::NumericMatrix& GXS,
    const Rcpp::NumericMatrix& GYS,
    const int n, const int p, const int q,
    const double lambda, const double psigx, const double psigy,
    const int n_iter, const double thresh,
    const Rcpp::CharacterVector& family)
{
    double * X = REAL(XS);
    double * Y = REAL(YS);
    double * GX = REAL(GXS);
    double * GY = REAL(GYS);
    std::string fam = Rcpp::as<std::string>(family);
    netreg::graph_penalized_linear_model_data data(X, Y,
                                                   GX, GY,
                                                   n, p, q,
                                                   lambda, 1.0,
                                                   psigx, psigy,
                                                   n_iter, thresh);

    netreg::edgenet* edge;
    if (fam == "binomial")
    {
        edge = new netreg::binomial_edgenet;
    }
    else if (fam == "gaussian")
    {
        edge = new netreg::gaussian_edgenet;
    }
    else
    {
        Rcpp::Rcerr << "No correct family given!" << "\n";
        return Rcpp::List::create(0);
    }
    edge->run(data);
    delete edge;
    return Rcpp::List::create(Rcpp::Named("coefficients") = data.coefficients(),
                              Rcpp::Named("intercept") = data.intercept());
}

