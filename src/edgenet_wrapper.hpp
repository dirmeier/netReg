/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_NETREGWRAPPER_HPP
#define NETREG_NETREGWRAPPER_HPP

/**
 * Wrapper for calling edge-net.
 *
 * Basically all the parameters are the same as in netReg.cpp gauss_edgenet
 */
void do_gauss_edgenet_(double *const X, double *const Y,
                       double *const GX, double *const GY,
                       int N, int P, int Q,
                       double lambda,
                       double psigyx, double psigy,
                       int n_iter, double thresh);

/**
 * Wrapper for calling cross-validation for edge-net.
 *
 * Basically all the parameters are the same as in netReg.cpp gauss_cv_edgenet
 */
void do_gauss_cv_edgenet_(double *const X, double *const Y,
                          double *const GX, double *const GY,
                          int n, int p, int q,
                          double psigyx, double psigy,
                          int n_iter, double thresh,
                          int n_folds, int *const foldid,
                          int n_foldid);

extern double *B_;
extern double *mu_;
extern double lamb_;
extern double psi_gx_;
extern double psi_gy_;
#endif //NETREG_NETREGWRAPPER_HPP
