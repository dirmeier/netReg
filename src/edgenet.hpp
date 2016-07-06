/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_EDGENET_H
#define NETREG_EDGENET_H

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "types.hpp"
#include "graph_penalized_linear_model_data.hpp"
#include "cv_fold.hpp"
#include "cv_set.hpp"

namespace netreg
{
    /**
     * Class for estimating the coeffiecients of a edge-regularized linear regression model.
     */
    class edgenet
    {
    public:
        /**
         * Calulates the coefficients of a graph-regularized regression model.
         *
         * @param modelData an object that holds all required data for the model
         */
        void run(graph_penalized_linear_model_data &modelData);

        matrix<double> mccd_
            (graph_penalized_linear_model_data &data,
             const double lambda,
             const double alpha,
             const double psigx,
             const double psigy,
             cv_fold &fold) const;

        /**
         * Calculate an univariate cyclic coordinate descent on the index of
         * a column of a multivariate response matrix.
         *
         * @param data an object containing all relevant data
         * @param B the current coefficient matrix
         * @param B_old the old coefficient matrix
         * @param qi the index of a column of the multivariate response amtrox
         */
        void uccd_
            (graph_penalized_linear_model_data &data,
             matrix<double> &B,
             matrix<double> &B_old,
             const int qi);

        void uccd_
            (graph_penalized_linear_model_data &data,
             matrix<double> &cfs,
             matrix<double> &o_cfs,
             const double lamb, const double alph,
             const double psigx, const double psigy,
             cv_fold &fold, int qi,
             matrix<double> &trainTXX,
             matrix<double> &trainTXY) const;

    private:
        /**
         * Updates the softthresholding parameter and the normalization
         * constant using the graph Laplacians.
         *
         * @param s the softthresholding parameter to be set
         * @param norm the normalization constant to be set
         * @param cfs the current estimate of the coefficients
         * @param LX the normalized laplacian of X
         * @param LY the normalized laplacian of Y
         * @param P the number of covariables
         * @param Q the number of responses
         * @param pi the current index of the column of X
         * @param qi the current index of the column of Y
         * @param psigx the penalty for the Laplacian of X
         * @param psigy the penalty for the Laplacien of Y
         */
        void graph_penalize
            (double &s, double &norm,
             const double psigx, const double psigy,
             matrix<double> &LX, matrix<double> &LY, matrix<double> &cfs,
             const int P, const int Q,
             const int pi, const int qi);

        /**
         * Adds the penalty from the Laplacian of X to the softthresholding
         * parameter s and the normalization constant norm.
         *
         * @param s the softthresholding parameter to be set
         * @param norm the normalization constant to be set
         * @param cfs the current estimate of the coefficients
         * @param LX the normalized laplacian of X
         * @param P the number of covariables
         * @param pi the current index of the column of X
         * @param qi the current index of the column of Y
         * @param psigx the penalty for the Laplacian of X
         */
        void lx_penalize
            (double &s, double &norm, const double psigx, matrix<double> &LX,
             matrix<double> &cfs, const int P, const int pi, const int qi);


        /**
         * Adds the penalty from the Laplacian of X to the softthresholding
         * parameter s and the normalization constant norm.
         *
         * @param s the softthresholding parameter to be set
         * @param norm the normalization constant to be set
         * @param psigy the penalty for the Laplacien of Y
         * @param LY the normalized laplacian of Y
         * @param B the current estimate of the coefficients
         * @param Q the number of responses
         * @param pi the current index of the column of X
         * @param qi the current index of the column of Y
         */
        void ly_penalize
            (double &s, double &norm, const double psigy,
             matrix<double> &LY, matrix<double> &B, const int Q,
             const int pi, const int qi);

        /**
         * Calculates the softhresholding parameter as well as the normalisation constant.
         *
         * @param s the softthresholding parameter to be set
         * @param norm the normalization constant to be set
         * @param TXX the square of the design matrix
         * @param TXY the design times the response matrix
         * @param B the current estimate of the coefficients
         * @param LX the normalized laplacian of X
         * @param LY the normalized laplacian of Y
         * @param P the number of covariables
         * @param Q the number of responses
         * @param pi the current index of the column of X
         * @param qi the current index of the column of Y
         * @param psigx the penalty for the Laplacian of X
         * @param psigy the penalty for the Laplacien of Y
         * @param lower boolean flag whether only the lower triangular matrix of TXX is initialized
         */
        void set_params
            (double &s, double &norm,
             matrix<double> &TXX, matrix<double> &TXY,
             matrix<double> &B,
             matrix<double> &LX, matrix<double> &LY,
             const int P, const int Q, const int pi, const int qi,
             const double psigx, const double psigy,
             const bool lower) const;

        /**
         * In cross-validation most of the parameters can be updated iteratively when needed.
         * This is done here and probably buggy, since a lot is calculated on the fly.
         * TODO check this
         *
         * @param soft reference to softthresholding parameter
         * @param norm reference to normalization term
         * @param TXX X'X matrix
         * @param TXY X'Y matrix
         * @param TXX_train X'X matrix only with training indeexs
         * @param TXY_train X'Y matrix only with testing indexes
         * @param train_idxs the indexes of the training set
         * @param test_idxs the indexes of the testing set
         * @param B the current estimate of the coefficients
         * @param X the design matrix
         * @param Y the response matrix
         * @param LX normalized graph laplacian for X
         * @param LY normalized graph laplacian for Y
         * @param pi current coefficient index
         * @param qi current response index
         * @param psigx penalization for LX
         * @param psigy penalization for LY
         * @param P number of covariables
         * @param Q number of responses
         */
        void preset_params
            (double &soft, double &norm,
             matrix<double> &TXX,
             matrix<double> &TXY,
             matrix<double> &TXX_train,
             matrix<double> &TXY_train,
             std::vector<int> &train_idxs, std::vector<int> &testIidxs,
             matrix<double> &B,
             matrix<double> &X, matrix<double> &Y,
             matrix<double> &LX, matrix<double> &LY,
             const int pi, const int qi,
             const double psigx, const double psigy,
             const int P, const int Q) const;

        /**
         * Update the normalization term incross-validation
         *
         * @param X design matrix
         * @param txx an element of X'X matrix
         * @param test_idxs the vector of testing indexes
         * @param cidx the current coefficient
         */
        double up_norm(matrix<double> &X,
                       const double txx,
                       std::vector<int> &test_idxs,
                       const int cidx);

        /**
         * Preset the TXX matrix at a specific element
         *
         * @param TXX X'X matrix
         * @param X the design matrix
         * @param cidx current coefficient index
         * @param row current row
         */
        void pre_txx(matrix<double> &TXX,
                     matrix<double> &X,
                     const int cidx,
                     const int row);
    };
}
#endif //NETREG_EDGENET_H
