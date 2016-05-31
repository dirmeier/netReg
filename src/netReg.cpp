/**
 * Author: Simon Dirmeier
 * Date: 03/04/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "netreg_wrapper.hpp"

/**
 * Extern C++ binding for R! (name mangling)
 */

double *B_;
double *mu_;

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
 * @param niters max number of iterations if parameter estimation does not converge in time
 * @param threshs convergence threshold
 */
SEXP gauss_edgenet(SEXP XS, SEXP YS,
             SEXP GXS, SEXP GYS,
             SEXP ns, SEXP ps, SEXP qs,
             SEXP lambdass,
             SEXP psi_gxs, SEXP psi_gys,
             SEXP niters, SEXP threshs)
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
    const double LAMBDA = (*REAL(lambdass));
    // cast R weighting values for GX to pointer to double
    const double PSI_GX = (*REAL(psi_gxs));
    // cast R weighting values for GY to pointer to double
    const double PSI_GY = (*REAL(psi_gys));
    // get number of max iterations
    const int N_ITER = (*INTEGER(niters));
    // get convergence threshold
    const double THRESH = (*REAL(threshs));

    // wrap the data
    run(X, Y, GX, GY,
                N, P, Q,
                LAMBDA, PSI_GX, PSI_GY,
                N_ITER,
                THRESH);
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
    // release SEXPs for garbage collection (what a stupid feature to need that, why not create vars on stack???)
    UNPROTECT(prtCnt);
    // return results to R
    return OS;
}

SEXP gauss_cv_edgenet(SEXP XS, SEXP YS, SEXP GXS, SEXP GYS,
                      SEXP ns, SEXP ps, SEXP qs,
                      SEXP lambs, SEXP psgxs, SEXP psgys,
                      SEXP niters, SEXP threshs,
                      SEXP nfolds, SEXP foldids)
{
    return R_NilValue;
}


};
