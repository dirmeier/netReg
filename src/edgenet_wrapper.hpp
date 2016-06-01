/**
 * Author: Simon Dirmeier
 * Date: 14/05/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_NETREGWRAPPER_HPP
#define NETREG_NETREGWRAPPER_HPP

void do_gauss_edgenet_(double *const X, double *const Y,
                       double *const GX, double *const GY,
                       const int N, const int P, const int Q,
                       const double LAMBDA, const double PSI_GX,
                       const double PSI_GY,
                       const int N_ITER, const double THRESH);

void do_gauss_cv_edgenet_(double *const X, double *const Y,
                          double *const GX, double *const GY,
                          const int N, const int P, const int Q,
                          const double LAMBDA, const double PSI_GX,
                          const double PSI_GY,
                          const int N_ITER, const double THRESH,
                          const int N_FOLDS, int * const foldid);

extern double *B_;
extern double *mu_;
extern double lamb_;
extern double psi_gx_;
extern double psi_gy_;
#endif //NETREG_NETREGWRAPPER_HPP
