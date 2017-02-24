/**
 * Author: Simon Dirmeier
 * Date: 24/02/17
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_EDGENET_GAUSSIAN_MODEL_SELECTION_HPP
#define NETREG_EDGENET_GAUSSIAN_MODEL_SELECTION_HPP

#include <map>
#include "graph_penalized_linear_model_cv_data.hpp"

namespace netreg
{

    class edgenet_gaussian_model_selection
    {
    public:
        /**
         * Find the set of optimal shrinkage parameters for a edge-penalized regression model.
         * Set is calculated using cross-validation.
         *
         * @param data
         *
         * @returns returns a map of shrinkage parameters
         */
        std::map<std::string, double> regularization_path
            (graph_penalized_linear_model_cv_data &data);
    };
}
#endif //NETREG_EDGENET_GAUSSIAN_MODEL_SELECTION_HPP
