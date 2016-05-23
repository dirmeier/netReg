/**
 * Author: Simon Dirmeier
 * Date: 14/05/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_GRAPH_FUNCTIONS_HPP
#define NETREG_GRAPH_FUNCTIONS_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "../struc/types.hpp"

namespace netreg
{
   /*
    * Calculate the normalized laplacian of a matrix.
    *
    * @param m the natrux fir which the laplacian is calculated
    * @return the normalized laplacian
    */
    matrix<double> laplacian(matrix<double> &m);

}
#endif //NETREG_GRAPH_FUNCTIONS_HPP
