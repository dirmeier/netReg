/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_PENALIZEDLINEARMODELDATA_HPP
#define NETREG_PENALIZEDLINEARMODELDATA_HPP

#include "linear_model_data.hpp"

#include <string>

namespace netreg
{
    /**
     * Data-structure for all required data for penalized regression.
     */
    class penalized_linear_model_data : public linear_model_data
    {

    public:
        /**
         * Constructor.
         *
         * @param x the nxp-dimensional design matrix of the linear model
         * @param y the nxq-dimensional response matrix of the linear model
         * @param n the number of samples (nrows X/Y)
         * @param p the number of covariables (ncols X)
         * @param q the number of responses (ncol Y)
         * @param niter max number of iterations in case estimation of the coefficients does not converge
         * @param thresh convergence threshold
         * @param is_column_first boolean if given matrices are column-first or not
         */
        penalized_linear_model_data(double *x, double *y,
                                    const int n, const int p, const int q,
                                    double const lambda, double const alpha,
                                    const int niter, const double thresh,
                                    std::string family)
            : linear_model_data(x, y, n, p, q, niter, thresh, family),
              ALPHA(alpha),
              LAMBDA(lambda)
        {
        }

        /**
         * Getter for the penalization parameter of the lasso.
         *
         * @return the penalization term for the LASSO
         */
        const double lambda()
        {
            return LAMBDA;
        }

        /**
         * Getter for the mixing parameter of the elastic-net.
         *
         * @return the mixnig term for the elastic net
         */
        const double alpha()
        {
            return ALPHA;
        }

    protected:
        const double ALPHA;   // mixing weights for elastic-net
        const double LAMBDA;  // penalization term for LASSO
    };
}
#endif //NETREG_PENALIZEDLINEARMODELDATA_HPP
