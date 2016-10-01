/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_EDGENETMODELSELECTION_HPP
#define NETREG_EDGENETMODELSELECTION_HPP

#include "graph_penalized_linear_model_cv_data.hpp"
#include "cv_set.hpp"

namespace netreg
{
    /**
     * Class that selects the best parameters of a model-selection using cross-validation
     */
    class edgenet_model_selection
    {
    public:

        /**
         * Calculate a regularization path for an edge-penalized regression model,
         * i.e. calculate different lambdas and psis.
         *
         * @param data the model data for which you want to estimate the
         *        optimal regularization parameters.
         *
         * @returns an pareto-optimal point
         */
        std::vector<double> regularization_path
            (graph_penalized_linear_model_cv_data &data);
    };
}
#endif //NETREG_EDGENETMODELSELECTION_HPP
