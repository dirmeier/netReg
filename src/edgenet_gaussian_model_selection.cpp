/**
 * Author: Simon Dirmeier
 * Date: 24/02/17
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include "edgenet_gaussian_model_selection.hpp"

#include "family.hpp"
#include "optim.hpp"
#include "edgenet_gaussian_loss_function.hpp"

namespace netreg
{

    std::map<std::string, double> edgenet_gaussian_model_selection::regularization_path(
        graph_penalized_linear_model_cv_data &data) const
    {
        optim opt;

        std::vector<double> start{0, 0, 0};
        std::vector<double> lower_bound{0.0, 0.0, 0.0};
        std::vector<double> upper_bound{100.0, 10000.0, 10000.0};
        const double rad_start = 0.49, rad_end = 1e-6;
        const int niter = 100000;

        switch (data.distribution_family())
        {
            case family::GAUSSIAN:
            default:
            {
                std::map<std::string, double> res =
                    opt.bobyqa<edgenet_gaussian_loss_function>(
                        data, start, lower_bound, upper_bound, rad_start,
                        rad_end, niter);
                return res;
            }
        }
    }
}
