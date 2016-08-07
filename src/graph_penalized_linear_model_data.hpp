/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_GRAPHPENALIZEDLINEARMODELDATA_HPP
#define NETREG_GRAPHPENALIZEDLINEARMODELDATA_HPP

#include "penalized_linear_model_data.hpp"

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "types.hpp"
#include "graph_functions.hpp"

namespace netreg
{
    /**
     * Data-structure for all required data for graph-penalized regression.
     */
    class graph_penalized_linear_model_data
        : public penalized_linear_model_data
    {

    public:

        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param gx a pxp-dimensional prior graph for the covariables of x
         * @param gy a qxq-dimensional prior graph for the responses of y
         * @param n the number of samples (nrows X/Y)
         * @param p the number of covariables (ncols X)
         * @param q the number of responses (ncol Y)
         * @param lambda a vector of length q of penalisation values for q univariate models
         * @param alpha a vector of length q of weightings for lasso/ridge
         * @param psi_gx a vector of length q of how much influence GX should have on the penalization
         * @param psi_gy a vector of length q of how much influence GY should have on the penalization
         * @param niter max number of iterations in case estimation of the coefficients does not converge
         * @param thresh convergence threshold
         */
        graph_penalized_linear_model_data(double *const x, double *const y,
                                          double *const gx, double *const gy,
                                          const int n, const int p,
                                          const int q,
                                          double const lambda,
                                          const double alpha,
                                          const double psi_gx,
                                          const double psi_gy,
                                          const int niter,
                                          const double thresh)
            : penalized_linear_model_data(x, y, n, p, q, lambda,
                                          alpha, niter, thresh),
              psi_gx(psi_gx), psi_gy(psi_gy),
              GX(gx, p, p), GY(gy, q, q),
              LX(laplacian(gx, p, p)), LY(laplacian(gy, q, q))
        {
        }

        /**
         * Getter for penalization term for laplacian of X
         *
         * @return the penalization term for X
         */
        const double psigx()
        {
            return psi_gx;
        }

        /**
         * Getter for penalization term for laplacian of Y
         *
         * @return the penalization term for Y
         */
        const double psigy()
        {
            return psi_gy;
        }

        /**
         * Getter for laplacian matrix of X
         *
         * @return reference to laplacian matrix
         */
        matrix<double> &lx()
        {
            return LX;
        }

        /**
         * Getter for laplacian matrix of X
         *
         * @return reference to laplacian matrix
         */
        matrix<double> &ly()
        {
            return LY;
        }
        const double psi_gx;  // Penalization vector for GX
        const double psi_gy;  // Penalization vector for GY
        matrix<double> GX;    // prior graph for the design matrix
        matrix<double> GY;    // prior graph for response matrix
        matrix<double> LX;    // Normalized Laplacian of GX
        matrix<double> LY;    // Normalized Laplacian of GY
    };
}
#endif //NETREG_GRAPHPENALIZEDLINEARMODELDATA_HPP
