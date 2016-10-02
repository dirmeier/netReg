/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include <memory>
#include <string>
#include <vector>

#include "graph_penalized_linear_model_data.hpp"
#include "graph_penalized_linear_model_cv_data.hpp"
#include "edgenet.hpp"
#include "edgenet_binomial.hpp"
#include "edgenet_gaussian.hpp"
#include "edgenet_model_selection.hpp"
#include "family.hpp"

#include "Rcpp.h"

extern "C"
{
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
 * @param niters max number of iterations if parameter estimation
 *        does not converge in time
 * @param threshs convergence threshold
 * @param fs family of distribution the response
 */
SEXP edgenet
    (SEXP XS, SEXP YS, SEXP GXS, SEXP GYS,
     SEXP lambda, SEXP psigx, SEXP psigy,
     SEXP niter, SEXP thresh,
     SEXP fs)
{
    BEGIN_RCPP;
    std::string fam = Rcpp::as<std::string>(fs);
    family f = fam == "binomial" ? family::BINOMIAL :
               fam == "gaussian" ? family::GAUSSIAN
                                 : family::NONE;
    if (f == family::NONE)
    {
        Rcpp::Rcerr << "Wrong family given!" << "\n";
        return R_NilValue;
    }
    Xdim = getAttrib(XS, R_DimSymbol);
    Ydim = getAttrib(YS, R_DimSymbol);
    netreg::graph_penalized_linear_model_data data
        (REAL(XS), REAL(YS), REAL(GXS), REAL(GYS),
         NTEGER(Xdim)[0], NTEGER(Xdim)[1], NTEGER(Ydim)[0],
         Rcpp::as<double>(lambda), 1.0,
         Rcpp::as<double>(psigx), Rcpp::as<double>(psigy),
         Rcpp::as<int>(niter), Rcpp::as<double>(thresh), f);
    // TODO change that back and include family in data
    netreg::edgenet edge;
    return edge.run(data);
    END_RCPP;
}
};
/**
 * Implementation of cross-validation for Edgenet.
 *
 * Finds and returns the optimal shrinkage values given a specific data-set.
 *
 * @param XS the (ns x ps)-dimensional design matrix
 * @param YS the (ns x qs)-dimensional response matrix
 * @param GXS the (ps x ps)-prior graph for XS
 * @param GYS the (qs x qs)-prior graph for YS
 * @param ns number of observations/samples/rows
 * @param ps number of covariates
 * @param qs number of responses
 * @param psi_gxs weighting value of GX
 * @param psi_gys weighting value of GY
 * @param niters max number of iterations if parameter estimation
 *        does not converge in time
 * @param threshs convergence threshold
 * @param nfolds the number of cross-validation sets created (as in k-fold cv)
 * @param foldids integer vector of assignments of observations
 *        to folds (i.e. vector of ns elements,  \in {1, ..., nfolds}
 * @param lenfoldids length of the vector above
 * @param familys family of distribution the response
 */
SEXP cv_edgenet(SEXP XS, SEXP YS, SEXP GXS, SEXP GYS,
                SEXP ns, SEXP ps, SEXP qs,
                SEXP psi_gxs, SEXP psi_gys,
                SEXP niters, SEXP threshs,
                SEXP nfolds, SEXP foldids, SEXP lenfoldids,
                SEXP familys)
{
    double *X = REAL(XS);
    // cast R response matrix to pointer to double
    double *Y = REAL(YS);
    // cast R prior matrix for X to pointer to double
    double *GX = REAL(GXS);
    // cast R prior matrix for Y to pointer to double
    double *GY = REAL(GYS);
    // get number of samples
    const int N = (*INTEGER(ns));
    // get number of covariables
    const int P = (*INTEGER(ps));
    // get number of responses
    const int Q = (*INTEGER(qs));
    // cast R design matrix to pointer to double
    // cast R weighting values for GX to pointer to double
    const double psigx = (*REAL(psi_gxs));
    // cast R weighting values for GY to pointer to double
    const double psigy = (*REAL(psi_gys));
    // get number of max iterations
    const int niter = (*INTEGER(niters));
    // get convergence threshold
    const double thresh = (*REAL(threshs));
    // the number of folds
    const int nfold = (*INTEGER(nfolds));
    // the fold assignments (usually not given)
    int *fold_ids = (INTEGER(foldids));
    // the length of the fold_ids array
    const int foldid_len = (*INTEGER(lenfoldids));
    // get family
    std::string family =
        CHAR(STRING_ELT(familys, 0)) == 'b' ? "binomial" : "gaussian";
    // call wrapper
    netreg::edgenet_model_selection e;
    std::vector<double> pop;
    // TODO fold parsing
    netreg::graph_penalized_linear_model_cv_data data;

    if (n_foldid == n)
        data = netreg::graph_penalized_linear_model_cv_data
            (X, Y, GX, GY, N, P, Q, -1, 1.0, psigx, psigy, niter, thresh,
             fold_ids, family);
    else
        data = netreg::graph_penalized_linear_model_cv_data
            (X, Y, GX, GY, N, P, Q, -1, 1.0, psigx, psigy, niter, thresh,
             nfold, family);
    std::vector<double> pop = e.regularization_path(data, fam);

    const double lamb_est = pop[0];
    const double psi_gx_est = psigx == -1 ? pop[1] : 0.0;
    const double psi_gy_est = psigy == -1 ? pop[2] : 0.0;
    // protection counter that is needed for not activating garbage collection
    int prtCnt = 0;
    SEXP shrink = PROTECT(allocVector(REALSXP, 3));
    prtCnt++;
    REAL(shrink)[0] = lamb_;
    REAL(shrink)[1] = psi_gx_;
    REAL(shrink)[2] = psi_gy_;
    // fold ids
    SEXP folds = PROTECT(allocVector(INTSXP, N));
    prtCnt++;
    std::vector &fold_ids = data.fold_ids();
    for (int i = 0; i < N; ++i) INTEGER(folds)[i] = fold_ids;
    // create a R list of size 2 that can be returned
    SEXP OS = PROTECT(allocVector(VECSXP, 2));
    prtCnt++;
    // set first element of list to the coef matrix
    SET_VECTOR_ELT(OS, 0, shrink);
    // set second element of list to intercept vector
    SET_VECTOR_ELT(OS, 1, folds);
    // create name array
    SEXP nms = PROTECT(allocVector(STRSXP, 2));
    prtCnt++;
    SET_STRING_ELT(nms, 0, mkChar("shrinkage_parameters"));
    SET_STRING_ELT(nms, 1, mkChar("fold_ids"));
    // assign names to list
    setAttrib(OS, R_NamesSymbol, nms);
    // release SEXPs for garbage collection
    UNPROTECT(prtCnt);
    // return results to R
    return OS;
}

