/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include <Rcpp.h>

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
// [[Rcpp
// [[Rcpp::export]]
Rcpp::List edgenet_(
    Rcpp::NumericMatrix XS, Rcpp::NumericMatrix YS,
    Rcpp::NumericMatrix GXS, Rcpp::NumericMatrix GYS,
    const int ns, const int ps, const int qs,
    Rcpp::NumericVector lambdass,
    Rcpp::NumericVector psi_gxs, Rcpp::NumericVector psi_gys,
    const int niters, const double threshs,
    Rcpp::CharacterVector familys)
{
    // get number of samples
    const int N = (*INTEGER(ns));
    // get number of covariables
    const int P = (*INTEGER(ps));
    // get number of responses
    const int Q = (*INTEGER(qs));
    // cast R design matrix to pointer to double
    double *X = REAL(XS);
    // cast R response matrix to pointer to double
    double *Y = REAL(YS);
    // cast R prior matrix for X to pointer to double
    double *GX = REAL(GXS);
    // cast R prior matrix for Y to pointer to double
    double *GY = REAL(GYS);
    // cast R lambda values to pointer to double
    const double lambda = (*REAL(lambdass));
    // cast R weighting values for GX to pointer to double
    const double psigx = (*REAL(psi_gxs));
    // cast R weighting values for GY to pointer to double
    const double psigy = (*REAL(psi_gys));
    // get number of max iterations
    const int niter = (*INTEGER(niters));
    // get convergence threshold
    const double thresh = (*REAL(threshs));
    // get family
    const char *family = CHAR(STRING_ELT(familys, 0));
    // call wrapper
    do_edgenet_(X, Y, GX, GY,
                N, P, Q,
                lambda, psigx, psigy,
                niter, thresh, family);
    // protection counter that is needed for not activating garbage collection
    int prtCnt = 0;
    // R object for coefficient matrix
    SEXP BS = PROTECT(allocMatrix(REALSXP, P, Q));
    // protect from gc
    prtCnt++;
    // R object for intercept vector
    SEXP intercept = PROTECT(allocVector(REALSXP, Q));
    // protect from gc
    prtCnt++;
    double *B = REAL(BS);
    double *b0 = REAL(intercept);
    for (int i = 0; i < Q; ++i)
    {
        // safe intercepts for R vector
        b0[i] = mu_[i];
        for (int j = 0; j < P; ++j)
        {
            // safe coefficients for R matrix
            B[j + P * i] = B_[j + P * i];
        }
    }
    // create a R list of size 2 that can be returned
    SEXP OS = PROTECT(allocVector(VECSXP, 2));
    prtCnt++;
    // set first element of list to the coef matrix
    SET_VECTOR_ELT(OS, 0, BS);
    // set second element of list to intercept vector
    SET_VECTOR_ELT(OS, 1, intercept);
    // create name array
    SEXP nms = PROTECT(allocVector(STRSXP, 2));
    prtCnt++;
    SET_STRING_ELT(nms, 0, mkChar("coefficients"));
    SET_STRING_ELT(nms, 1, mkChar("intercept"));
    // assign names to list
    setAttrib(OS, R_NamesSymbol, nms);
    // release SEXPs for garbage collection
    UNPROTECT(prtCnt);
    // return results to R
    return OS;
}

