/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "edgenet.hpp"
#include "math_functions.hpp"
#include "stat_functions.hpp"

namespace netreg
{
    void edgenet::run(graph_penalized_linear_model_data &data)
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();
        matrix<double> &coef = data.coefficients();
        matrix<double> old_coef(P, Q);
        cvector<double> &intr = data.intercept();
        const double THRESH = data.threshold();
        const int N_ITER = data.max_iter();
        int iter = 0;
        do
        {
            // TODO is this parallizable, guess not
            for (int qi = 0; qi < Q; ++qi)
                uccd_(data, coef, old_coef, qi);
        }

        while (arma::accu(arma::abs(coef - old_coef)) > THRESH &&
               iter++ < N_ITER);
        // calculate intercepts of the linear model
        intr = intercept(data.design(), data.response(), coef);
    }

    void edgenet::uccd_
        (graph_penalized_linear_model_data &data,
         matrix<double> &coef,
         matrix<double> &old_coef,
         const int qi)
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();
        const double THRESH = data.threshold();
        const int N_ITER = data.max_iter();
        const double LAMBDA = data.lambda();
        const double ALPHA = data.alpha();
        const double PSI_GX = data.psigx();
        const double PSI_GY = data.psigy();
        matrix<double> &TXX = data.txx();
        matrix<double> &TXY = data.txy();
        matrix<double> &LX = data.lx();
        matrix<double> &LY = data.ly();
        // weighted penalization param of Elastic-net
        const double lalph = ALPHA * LAMBDA;
        // normalization for soft-thresholding
        const double enorm = 1.0 + LAMBDA * (1 - ALPHA);
        // iteration counter
        int iter = 0;
        // do while estimation of params does not converge
        do
        {
            // fix bnew_i and calculate least-squares
            // coefficient on partial residual
            for (int pi = 0; pi < P; ++pi)
            {
                // safe current estimate of coefficients
                old_coef(pi, qi) = coef(pi, qi);
                double s = 0.0;
                double norm = 0.0;
                set_params
                    (s, norm, TXX, TXY, coef,
                     LX, LY, P, Q, pi, qi, PSI_GX, PSI_GY, false);
                // soft-thresholded version of estimate
                coef(pi, qi) = softnorm(s, lalph, enorm * norm);
            }
        }
        while (
            arma::accu(arma::abs(coef.col(qi) - old_coef.col(qi))) > THRESH &&
            iter++ < N_ITER);
    }

    void edgenet::uccd_
        (const int P, const int Q,
         const double thresh, const int niter,
         const double lambda, const double alpha,
         const double psigx, const double psigy,
         matrix<double> &TXX, matrix<double> &TXY,
         matrix<double> &LX, matrix<double> &LY,
         matrix<double> &coef,
         matrix<double> &old_coef,
         const int qi) const
    {


//        std::cout << std::endl;
//        for (int i = 0; i < TXX.n_rows; ++i)
//        {
//            for (int j = 0; j < TXX.n_cols; ++j)
//            {
//                std::cout << TXX(i, j) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//        for (int i = 0; i < TXY.n_rows; ++i)
//        {
//            for (int j = 0; j < TXY.n_cols; ++j)
//            {
//                std::cout << TXY(i, j) << " ";
//            }
//            std::cout <<std::endl;
//        }
//        std::cout <<  std::endl;

        // weighted penalization param of Elastic-net
        const double lalph = alpha * lambda;
        // normalization for soft-thresholding
        const double enorm = 1.0 + lambda * (1 - alpha);
        // iteration counter
        int iter = 0;
        // do while estimation of params does not converge
        do
        {
            // fix bnew_i and calculate least-squares
            // coefficient on partial residual
            for (int pi = 0; pi < P; ++pi)
            {
                // safe current estimate of coefficients
                old_coef(pi, qi) = coef(pi, qi);
                double s = 0.0;
                double norm = 0.0;
                set_params
                    (s, norm, TXX, TXY, coef,
                     LX, LY, P, Q, pi, qi, psigx, psigy, false);
                // soft-thresholded version of estimate
                coef(pi, qi) = softnorm(s, lalph, enorm * norm);
            }
        }
        while (
            arma::accu(arma::abs(coef.col(qi) - old_coef.col(qi))) > thresh &&
            iter++ < niter);
    }

    void edgenet::set_params
        (double &s, double &norm,
         matrix<double> &TXX, matrix<double> &TXY,
         matrix<double> &coef,
         matrix<double> &LX,
         matrix<double> &LY,
         const int P, const int Q,
         const int pi, const int qi,
         const double psigx,
         const double psigy,
         const bool lower) const
    {
        s = pls(TXX, TXY, coef, pi, qi, P, lower);
        norm = (TXX)(pi, pi);
        graph_penalize(s, norm, psigx, psigy,
                       LX, LY, coef,
                       P, Q, pi, qi);
    }

    void edgenet::graph_penalize
        (double &s, double &norm,
         const double psigx, const double psigy,
         matrix<double> &LX, matrix<double> &LY, matrix<double> &cfs,
         const int P, const int Q,
         const int pi, const int qi) const
    {
        if (psigx != 0)
            lx_penalize(s, norm, psigx, LX, cfs, P, pi, qi);
        if (psigy != 0)
            ly_penalize(s, norm, psigy, LY, cfs, Q, pi, qi);
    }

    void edgenet::lx_penalize
        (double &s, double &norm, const double psigx,
         matrix<double> &LX, matrix<double> &cfs, const int P,
         const int pi, const int qi) const
    {
        if (psigx == 0)
            return;
        double xPenalty = -LX(pi, pi) * cfs(pi, qi);
        for (int j = 0; j < P; j++)
            xPenalty += LX(pi, j) * cfs(j, qi);
        s = s - 2 * psigx * xPenalty;
        norm += 2 * psigx * LX(pi, pi);
    }

    void edgenet::ly_penalize
        (double &s, double &norm, const double psigy,
         matrix<double> &LY, matrix<double> &cfs, const int Q,
         const int pi, const int qi) const
    {
        if (psigy == 0 || Q == 1)
            return;
        // penalization for GY
        double yPenalty = -cfs(pi, qi) * LY(qi, qi);
        for (int j = 0; j < Q; ++j)
            yPenalty += cfs(pi, j) * LY(j, qi);
        s = s - 2 * psigy * yPenalty;
        norm += 2 * psigy * LY(qi, qi);
    }

    matrix<double> edgenet::mccd_(graph_penalized_linear_model_data &data,
                                  const double lambda, const double alpha,
                                  const double psigx, const double psigy,
                                  cv_fold &fold) const
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();
        matrix<double> coef(P, Q, arma::fill::ones);
        matrix<double> old_coef(P, Q);
        const double thresh = data.threshold();
        const int niter = data.max_iter();
        matrix<double> &X = data.design();
        matrix<double> &Y = data.response();
        matrix<double> &LX = data.lx();
        matrix<double> &LY = data.ly();
        index_vector &trainIdxs = fold.train_set();
        index_vector &testIdxs = fold.test_set();
        matrix<double> Xtrain = X.rows(trainIdxs);
        matrix<double> Ytrain = Y.rows(trainIdxs);
        matrix<double> TXtrain = Xtrain.t();
        matrix<double> train_txx = TXtrain * Xtrain;
        matrix<double> train_txy = TXtrain * Ytrain;

//        std::cout <<  "full" << std::endl;
//        for (int i = 0; i < data.txx().n_rows; ++i)
//        {
//            for (int j = 0; j < data.txx().n_cols; ++j)
//            {
//                std::cout << data.txx()(i, j) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;
//        for (int i = 0; i < Ytrain.n_rows; ++i)
//        {
//            for (int j = 0; j < Ytrain.n_cols; ++j)
//            {
//                std::cout << Ytrain(i, j) << " ";
//            }
//            std::cout << std::endl;
//        }
//        std::cout << std::endl;

        int iter = 0;
        do
        {
            // TODO is this parallizable? think not
            for (int qi = 0; qi < Q; ++qi)
            {
                uccd_(P, Q, thresh, niter, lambda, alpha, psigx, psigy,
                      train_txx, train_txy, LX, LY, coef, old_coef, qi);
            }
        }
        while (arma::accu(arma::abs(coef - old_coef)) > thresh &&
               iter++ < niter);

        return coef;
    }

    void edgenet::uccd_
        (graph_penalized_linear_model_data &data,
         matrix<double> &coef, matrix<double> &old_coef,
         const double lamb, const double alph, const double psigx,
         const double psigy, cv_fold &fold, int qi,
         matrix<double> &trainTXX, matrix<double> &trainTXY) const
    {
        const int P = data.covariable_count();
        const int Q = data.response_count();
        const double THRESH = data.threshold();
        const int N_ITER = data.max_iter();
        matrix<double> &X = data.design();
        matrix<double> &Y = data.response();
        matrix<double> &TXX = data.txx();
        matrix<double> &TXY = data.txy();
        matrix<double> &LX = data.lx();
        matrix<double> &LY = data.ly();
        // weighted penalization param of Elastic-net
        const double lalph = alph * lamb;
        // normalization for soft-thresholding
        const double enorm = 1.0 + lamb * (1 - alph);
        // iteration counter
        // the vector of training indexes of the current fold
        index_vector &trainIdxs = fold.train_set();
        // the vector of testing indexed of the current fold
        index_vector &testIdxs = fold.test_set();

        matrix<double> train_txx;
        matrix<double> train_txy;

        // boolean whether we calculate everything by hand or take the precomputed value
        bool isComputed = false;
        int iter = 0;
        // do while estimation of params does not converge
        do
        {
            // fix bnew_i and calculate least-squares
            // coefficient on partial residual
            for (int pi = 0; pi < P; ++pi)
            {
                old_coef(pi, qi) = coef(pi, qi);
                double norm = 0.0;
                double s = 0.0;
                if (!isComputed)
                {
//                    // make some docu here
////                    // computes s and norm and sets up matrices
//                    preset_params(s, norm, TXX, TXY, trainTXX, trainTXY,
//                                  trainIdxs, testIdxs, coef, X, Y, LX, LY,
//                                  pi, qi, psigx, psigy, P, Q);
//                    isComputed = pi == P - 1;
//                    std::cout << "\n";
                    for (int i = 0; i < P; ++i)
                    {
                        for (int j = 0; j < Q; ++j)
                        {
                            std::cout << trainTXY(i, j) << " ";
                        }
                        std::cout << "\n";
                    }
                    std::cout << "\n";

                }
//                else
//                    set_params
//                        (s, norm, trainTXX, trainTXY, coef,
//                         LX, LY, P, Q, pi, qi, psigx, psigy, true);

                // soft-thresholded version of estimate
                coef(pi, qi) = softnorm(s, lalph, enorm * norm);
            }
        }
        while (
            arma::accu(arma::abs(coef.col(qi) - old_coef.col(qi))) > THRESH &&
            iter++ < N_ITER);
    }

    void edgenet::preset_params
        (double &soft, double &norm,
         matrix<double> &TXX, matrix<double> &TXY,
         matrix<double> &trainTXX, matrix<double> &trainTXY,
         std::vector<int> &trainIdxs, std::vector<int> &testIdxs,
         matrix<double> &cfs, matrix<double> &X, matrix<double> &Y,
         matrix<double> &LX, matrix<double> &LY,
         const int pi, const int qi,
         const double psigx, const double psigy,
         const int P, const int Q) const
    {
        /*
        * Caution!
        *
        * This calculates 'soft', 'norm' and fills the matrices TXX and TXYqi.
        * However TXX is only completel filled after P iterations!!! So this step is error-prone.
        */
        // Iterate over all training indexes.
        for (auto &row : trainIdxs)
        {
            // Compute dot product of X'Y[cidx,qi]
            // // This matrix<double>is filled after one full cycle!!!!
            trainTXY(pi, qi) += X(row, pi) * Y(row, qi);
            // compute X'X[cidx, didx]
            // Only compute the diagonal and the lower half for speed reasons.
            // Caution: ALWAYS take trainTXX[max(i,j), min(i,j)]
            // -> possible cause X'X is symmetric
            // This matrix<double>is filled after one full cycle!!!!
            pre_txx(trainTXX, X, pi, row);
            // soft-thresholded value for bnew_i
            soft += X(row, pi) *
                    partial_residual(X, Y, cfs, row, pi, qi);
        }
        // Safe, because neither trainTXX nor trainTXYqi are required
        norm = up_norm(X, TXX(pi, pi), testIdxs, pi);
        // set up normalizer and update s
        graph_penalize(soft, norm, psigx, psigy, LX, LY, cfs, P, Q, pi, qi);
    }

    double edgenet::up_norm(matrix<double> &X, const double txx,
                            std::vector<int> &testIdxs,
                            const int cidx) const
    {
        double norm = txx;
        if (!testIdxs.size())
            return txx;
        for (auto &row : testIdxs)
            norm -= X(row, cidx) * X(row, cidx);
        return norm;
    }

    void edgenet::pre_txx(matrix<double> &txx, matrix<double> &X,
                          const int cidx, const int row) const
    {
        for (int didx = 0; didx <= cidx; ++didx)
            txx(cidx, didx) += X(row, cidx) * X(row, didx);
    }

}