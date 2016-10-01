/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "edgenet_model_selection.hpp"

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include <numeric>
#include <vector>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "edgenet_gaussian_loss_function.hpp"
#include "edgenet_binomial_loss_function.hpp"
#include "optim.hpp"

namespace netreg
{
    std::vector<double>
    edgenet_model_selection::regularization_path
        (graph_penalized_linear_model_data &data, const int nfolds,
         const char *fam)
    {
        cv_set cvset(data.design().n_rows, nfolds);
        return regularization_path_(data, cvset, folds, fam);
    }

    std::vector<double>
    edgenet_model_selection::regularization_path
        (graph_penalized_linear_model_data &data, int *const foldid,
         const char *fam)
    {
        cv_set cvset(data.design().n_rows, foldid);
        return regularization_path_(data, cvset, folds, fam);
    }

    std::vector<double>
    edgenet_model_selection::regularization_path_
        (graph_penalized_linear_model_data &data, cv_set &cvset,
         const char *fam)
    {
        set_fold_ids(data, cvset);
        optim opt;
        std::vector<double> start{0, 0, 0};
        std::vector<double> lower_bound{0.0, 0.0, 0.0};
        std::vector<double> upper_bound{100.0, 10000.0, 10000.0};
        const double rad_start = 0.49, rad_end = 1e-6;
        const int niter = 1000;
        switch (fam[0])
        {
            case 'b':
                return opt.bobyqa<binomial_edgenet_loss_function>
                              (data, cvset,
                               start, lower_bound, upper_bound,
                               rad_start, rad_end, niter);
            case 'g':
            default:
                return opt.bobyqa<gaussian_edgenet_loss_function>
                              (data, cvset,
                               start, lower_bound, upper_bound,
                               rad_start, rad_end, niter);

        }

    }

    void edgenet_model_selection::set_foldids
        (graph_penalized_linear_model_data &data, cv_set &cvset)
    {
        #pragma omp parallel for
        for (int i = 0; i < cvset.fold_count(); i++)
        {
            cv_fold &fold = cvset.get_fold(i);
            for (arma::uvec::iterator j = fold.test_set().begin();
                 j != fold.test_set().end(); ++j)
                folds[*j] = i;
        }
    }
}