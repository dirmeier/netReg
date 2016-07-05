/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */
#ifndef NETREG_OPTIM_HPP
#define NETREG_OPTIM_HPP

#include <string>
#include <vector>

#include "../inst/dlib/optimization.h"
#include "pareto_optimal_point.hpp"
#include "graph_penalized_linear_model_data.hpp"
#include "cv_set.hpp"
#include "edgenet_loss_function.hpp"

namespace netreg
{
    /**
     * Class that uses several optimization functions from the dlib library.
     * dlib <3
     */
    class optim
    {
    public:
        /**
         * Calculate the pareto optimal points for a given loss function using the BOBYQA algorithm by Powell
         * (MJD Powell, The BOBYQA algorithm for bound constrained optimization without derivatives, 2009)
         *
         * BOBYQA allows for box constraints and does not need derivates.
         * It finds a local optimum for non-convex loss functions defined by a specified radius.
         *
         * We use it to minimize functions from the loss_functions package
         *
         * @template loss_function the class-name of an objective function that should be minimized
         * @param point an intial parameter setting for bobyqa
         * @param data the model data for the loss function
         * @param cvset the cv-set on which the estimated predictor is tested against
         * @param lower_bound lower bound box constraint
         * @param upper_bound upper bound box constraint
         * @param radius_start initial value of trust region
         * @param radius_start final value of trust region
         * @param niter maximum calls to the loss function
         */
        template<typename loss_function>
        void bobyqa
            (pareto_optimal_point<std::string, double> &point,
             graph_penalized_linear_model_data &data, cv_set &cvset,
             std::vector<double> &lower_bound,
             std::vector<double> &upper_bound,
             const double radius_start, const double radius_stop,
             const int niter)
        {
            int sz = static_cast<int>(point.npar());
            // convert to dlib objects
            dlib::matrix<double> par(sz, 1), lb(sz, 1), ub(sz, 1);
            for (int i = 0; i < sz; ++i)
            {
                par(i, 0) = point[i].second;
                lb(i, 0) = lower_bound[i];
                ub(i, 0) = upper_bound[i];
            }
            // minimize the loss_function
            dlib::find_min_bobyqa(
                loss_function(data, cvset),
                par,
                par.size() * 2 + 1,
                lb,
                ub,
                radius_start,
                radius_stop,
                niter
            );
            for (int i = 0; i < sz; ++i)
            {
                point[i].second = par(i, 0);
            }
        }
    };
}
#endif //NETREG_OPTIM_HPP
