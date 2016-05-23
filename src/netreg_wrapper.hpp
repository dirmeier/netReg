/**
 * Author: Simon Dirmeier
 * Date: 14/05/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_NETREGWRAPPER_HPP
#define NETREG_NETREGWRAPPER_HPP

void run(double *const X, double *const Y,
         double *const GX, double *const GY,
         const int N, const int P, const int Q,
         const double LAMBDA, const double PSI_GX, const double PSI_GY,
         const int N_ITER, const double THRESH);

extern double *B_;
extern double *mu_;

#endif //NETREG_NETREGWRAPPER_HPP
