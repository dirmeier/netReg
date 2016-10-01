/**
 * Author: Simon Dirmeier
 * Date: 02/10/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include "edgenet.hpp"

#include "edgenet_binomial.hpp"
#include "edgenet_gaussian.hpp"

#include <string>

namespace netreg
{
    virtual void edgenet::run(graph_penalized_linear_model_data &data) const
    {
        if (data.family() == "binomial")
        {
            netreg::edgenet_binomial edge;
            edge.run(data);
        }
        else
        {
            netreg::edgenet_gaussian edge;
            edge.run(data);
        }
    }
}