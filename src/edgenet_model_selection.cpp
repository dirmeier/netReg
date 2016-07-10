/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "edgenet_model_selection.hpp"

#include <numeric>
#include <vector>

#include "edgenet_loss_function.hpp"
#include "optim.hpp"

namespace netreg
{
    pareto_optimal_point<std::string, double> edgenet_model_selection::regularization_path
        (graph_penalized_linear_model_data &data,
         const int nfolds)
    {
        cv_set cvset(data.design().n_rows, nfolds);
        return regularization_path(data, cvset);
    }

    pareto_optimal_point<std::string, double> edgenet_model_selection::regularization_path
        (graph_penalized_linear_model_data &data, cv_set &cvset)
    {
        optim opt;
        pareto_optimal_point<std::string, double> pop;
        pop.put("lambda", 100.0);;
        pop.put("psigx", 1000.0);
        pop.put("psigy", 1000.0);
        std::vector<double> lower_bound{0.0, 0.0, 0.0};
        std::vector<double> upper_bound{100.0, 1000.0, 1000.0};
        const double rad_start = 0.49, rad_end = 1e-6;
        const int niter = 1000;
        opt.bobyqa<edgenet_loss_function>(pop, data, cvset,
                                          lower_bound, upper_bound,
                                          rad_start, rad_end, niter);
        return pop;
    }
}