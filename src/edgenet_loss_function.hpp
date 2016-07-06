#ifndef NETREG_EDGENETLOSSFUNCTION_HPP
#define NETREG_EDGENETLOSSFUNCTION_HPP

#include <numeric>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

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
    class edgenet_loss_function
    {
    public:
        /**
         * Creates an objective function object that can be used for minimization using dlib.
         *
         * @param data the complete dataset required for edge-regularized regression
         * @param cvset a cross-validation set
         */
        edgenet_loss_function
            (graph_penalized_linear_model_data &data, cv_set &cvset) :
            data_(data),
            cvset_(cvset),
            X_(data_.design()),
            Y_(data.response()),
            nfolds_(static_cast<unsigned long>(cvset.fold_count()))
        { }

        /**
         * Over-write operator () in order to get functor functionality (object behaves like a function)
         *
         * @param params the free parameters on the objective function
         */
        double operator()(const dlib::matrix<double> &params) const
        {
            std::vector<double> sses(nfolds_);
            // do n-fold cross-validation
            for (int fc = 0; fc < nfolds_; ++fc)
            {
                cv_fold &fold = cvset_.get_fold(fc);
                matrix<double> coef =
                    edgenet_loss_function::edgenet_.mccd_(data_,
                                                          params(0),
                                                          params(1),
                                                          params(2),
                                                          params(3),
                                                          fold);
                double err = sse(coef, X_, Y_, fold.test_set());
                sses[fc] = err;
            }
            return std::accumulate(sses.begin(), sses.end(), 0.0)
                   / sses.size();
        }

    private:
        // data required for a edge-regularized regression model
        graph_penalized_linear_model_data &data_;
        cv_set &cvset_;          // cv-set on which the selected model is evaluated
        matrix<double> &X_;      // design matrix
        matrix<double> &Y_;      // response matrix
        int nfolds_;             // number of folds
        static constexpr edgenet edgenet_ = edgenet(); // static object so that re-instatiation is avoided
    };
}
#endif //NETREG_EDGENETLOSSFUNCTION_HPP
