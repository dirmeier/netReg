#ifndef NETREG_EDGENETLOSSFUNCTION_HPP
#define NETREG_EDGENETLOSSFUNCTION_HPP

#include <numeric>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#ifndef NDEBUG
#include <iostream>
#else
#include <omp.h>
#endif
#include "types.hpp"
#include "edgenet.hpp"
#include "graph_penalized_linear_model_data.hpp"
#include "cv_set.hpp"
#include "error_functions.hpp"
#include "../inst/dlib/matrix.h"

namespace netreg
{
    /**
     * Functor class representing the objective function of a edge-regularized regression model.
     */
    class gaussian_edgenet_loss_function
    {
    public:
        /**
         * Creates an objective function object that can be used for minimization using dlib.
         *
         * @param data the complete dataset required for edge-regularized regression
         * @param cvset a cross-validation set
         */
        gaussian_edgenet_loss_function
            (graph_penalized_linear_model_data &data,
             cv_set &cvset):
            data_(data),
            cvset_(cvset),
            X_(data.design()),
            Y_(data.response()),
            nfolds_(static_cast<int>(cvset.fold_count())),
            edgenet_(),
            do_psigx_(data.psigx() == -1),
            do_psigy_(data.psigy() == -1)
        { }

        /**
         * Over-write operator () in order to get functor functionality (object behaves like a function)
         *
         * @param params the free parameters on the objective function
         */
        double operator()(const dlib::matrix<double> &p) const
        {
            std::vector<double> sses(nfolds_);
            // do n-fold cross-validation
#ifdef NDEBUG
            #pragma omp parallel for
#endif
            for (int fc = 0; fc < nfolds_; ++fc)
            {
                cv_fold &fold = cvset_.get_fold(fc);
                matrix<double> coef;
                if (do_psigx_ && do_psigy_)
                {
                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, p(1, 0),
                                          p(2, 0), fold);
                }
                else if (do_psigy_)
                {
                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, 0, p(2, 0),
                                          fold);
                }
                else if (do_psigx_)
                {
                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, p(1, 0), 0,
                                          fold);
                }
                else
                {
                    coef = edgenet_.run_cv(data_, p(0, 0), 1.0, 0, 0, fold);
                }
#ifndef NDEBUG
                std::cout << " test idx ";
                for(arma::uvec::iterator i = fold.test_set().begin(); i != fold.test_set().end(); ++i){
                    std::cout << *i <<  " ";
                }
                std::cout << std::endl;
                std::cout <<  "coef " <<std::endl;
                for (int i = 0; i < coef.n_rows; i++){
                    for (int j = 0; j < coef.n_cols; j++){
                        std::cout << coef(i,j) << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
                std::cout << std::endl;
#endif
                double err = sse(coef, X_, Y_, fold.test_set());
#ifndef NDEBUG
                std::cout << "lambda " <<  p(0, 0) << " err " << err << std::endl;
#endif
                sses[fc] = err;
            }
            double  err = std::accumulate(sses.begin(), sses.end(), 0.0);
#ifndef NDEBUG
            std::cout << "lambda " <<  p(0, 0) << " final err " << err << std::endl << std::endl;
#endif
            return err;
        }

    private:
        // data required for a edge-regularized regression model
        graph_penalized_linear_model_data &data_;
        cv_set &cvset_;          // cv-set on which the selected model is evaluated
        matrix<double> &X_;      // design matrix
        matrix<double> &Y_;      // response matrix
        int nfolds_;             // number of folds
        const edgenet edgenet_;
        const bool do_psigx_;
        const bool do_psigy_;
    };
}
#endif //NETREG_EDGENETLOSSFUNCTION_HPP
