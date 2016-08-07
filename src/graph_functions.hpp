/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_GRAPH_FUNCTIONS_HPP
#define NETREG_GRAPH_FUNCTIONS_HPP

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "types.hpp"

namespace netreg
{
   /*
    * Calculate the normalized laplacian of a matrix.
    *
    * @param x the pointer for which the laplacian is calculated (col first)
    * @param n nrows of x
    * @param m ncols of y
    * @return the normalized laplacian
    */
    matrix<double> laplacian(const double * x, int n, int m);

}
#endif //NETREG_GRAPH_FUNCTIONS_HPP
