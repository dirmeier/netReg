/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "edgenet_model_selection.hpp"

#include <numeric>
#include <vector>
#include <iostream>


#include "edgenet_loss_function.hpp"
#include "optim.hpp"

namespace netreg
{
    std::vector<double>
    edgenet_model_selection::regularization_path
        (graph_penalized_linear_model_data &data, const int nfolds)
    {
        cv_set cvset(data.design().n_rows, nfolds);
        return regularization_path_(data, cvset);
    }

    std::vector<double>
    edgenet_model_selection::regularization_path
        (graph_penalized_linear_model_data &data, int *const foldid)
    {
        cv_set cvset(data.design().n_rows, foldid);
        return regularization_path_(data, cvset);
    }

    std::vector<double>
    edgenet_model_selection::regularization_path_
        (graph_penalized_linear_model_data &data, cv_set &cvset)
    {
        optim opt;
        std::vector<double> start{100.0, 10000.0, 10000.0};
        std::vector<double> lower_bound{0.0, 0.0, 0.0};
        std::vector<double> upper_bound{100.0, 10000.0, 10000.0};
        const double rad_start = 0.49, rad_end = 1e-6;
        const int niter = 1000;
        return opt.bobyqa<edgenet_loss_function>
                      (data, cvset,
                       start, lower_bound, upper_bound,
                       rad_start, rad_end, niter);
    }
}