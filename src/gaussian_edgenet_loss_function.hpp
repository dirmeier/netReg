#ifndef NETREG_EDGENETLOSSFUNCTION_HPP
#define NETREG_EDGENETLOSSFUNCTION_HPP

#include <numeric>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>
#include <omp.h>
#include "types.hpp"
#include "gaussian_edgenet.hpp"
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
            #pragma omp parallel for
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
                double err = sse(coef, X_, Y_, fold.test_set());
                sses[fc] = err;
            }
            double  err = std::accumulate(sses.begin(), sses.end(), 0.0);
            return err;
        }

    private:
        // data required for a edge-regularized regression model
        graph_penalized_linear_model_data &data_;
        cv_set &cvset_;          // cv-set on which the selected model is evaluated
        matrix<double> &X_;      // design matrix
        matrix<double> &Y_;      // response matrix
        int nfolds_;             // number of folds
        const gaussian_edgenet edgenet_;
        const bool do_psigx_;
        const bool do_psigy_;
    };
}
#endif //NETREG_EDGENETLOSSFUNCTION_HPP
