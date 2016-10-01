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
        (graph_penalized_linear_model_cv_data &data)
    {
        optim opt;
        std::vector<double> start{0, 0, 0};
        std::vector<double> lower_bound{0.0, 0.0, 0.0};
        std::vector<double> upper_bound{100.0, 10000.0, 10000.0};
        const double rad_start = 0.49, rad_end = 1e-6;
        const int niter = 1000;
        switch (data.family())
        {
            case 'b':
                return opt.bobyqa<edgenet_binomial_loss_function>
                              (data, start, lower_bound, upper_bound,
                               rad_start, rad_end, niter);
            case 'g':
            default:
                return opt.bobyqa<edgenet_gaussian_loss_function>
                              (data, start, lower_bound, upper_bound,
                               rad_start, rad_end, niter);

        }
    }
}