/**
 * Author: Simon Dirmeier
 * Date: 01/10/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include "graph_penalized_linear_model_cv_data.hpp"

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <vector>
#include <cmath>

namespace netreg
{
    void graph_penalized_linear_model_cv_data::set_fold_ids()
    {
        if (static_cast<int>(fold_ids_.size()) != N)
            fold_ids_.resize(static_cast<std::vector<int>::size_type>(N));
        #pragma omp parallel for
        for (int i = 0; i < cvset_.fold_count(); i++)
        {
            cv_fold &fold = cvset_.get_fold(i);
            for (arma::uvec::iterator j = fold.test_set().begin();
                 j != fold.test_set().end(); ++j)
                fold_ids_[*j] = i;
        }
    }
}