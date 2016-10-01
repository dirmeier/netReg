/**
 * Author: Simon Dirmeier
 * Date: 01/10/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_GRAPH_PENALIZED_LINEAR_MODEL_CV_DATA_HPP
#define NETREG_GRAPH_PENALIZED_LINEAR_MODEL_CV_DATA_HPP

#include "graph_penalized_linear_model_data.hpp"

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "types.hpp"

namespace netreg
{
    /**
     * Data-structure for all required data for graph-penalized regression.
     */
    class graph_penalized_linear_model_cv_data
        : public graph_penalized_linear_model_data
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
         * @param nfolds the number of folds
         */
        graph_penalized_linear_model_cv_data
            (double *const x, double *const y,
             double *const gx, double *const gy,
             const int n, const int p, const int q,
             double const lambda, const double alpha,
             const double psi_gx, const double psi_gy,
             const int niter, const double thresh,
             const int nfolds)
            : graph_penalized_linear_model_data
                  (x, y, gx, gy, n, p, q, lambda, alpha, niter, thresh),
              fold_ids_(design().n_rows)
        {
            cvset_ = cv_set(design().n_rows, fold_ids);
            set_fold_ids();
        }

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
         * @param fold_ids fold id mappings
         */
        graph_penalized_linear_model_cv_data
            (double *const x, double *const y,
             double *const gx, double *const gy,
             const int n, const int p, const int q,
             double const lambda, const double alpha,
             const double psi_gx, const double psi_gy,
             const int niter, const double thresh,
             const int* fold_ids)
            : graph_penalized_linear_model_data
                  (x, y, gx, gy, n, p, q, lambda, alpha, niter, thresh),
              fold_ids_(design().n_rows)
        {
            cvset_ = cv_set(design().n_rows, fold_ids);
            set_fold_ids();
        }

    private:
        void set_fold_ids();
        // mapping from fold id to index in samples
        std::vector<int> fold_ids_;
        // the cross validation folds
        cv_set cvset_;
    };
}

#endif //NETREG_GRAPH_PENALIZED_LINEAR_MODEL_CV_DATA_HPP
