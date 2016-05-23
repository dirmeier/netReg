#ifndef NETREG_TYPES_HPP
#define NETREG_TYPES_HPP

#include <armadillo>

/**
 * Some typedefs for easier readability/writability.
 */
template<class T>
using matrix = arma::Mat<T>;
template<class T>
using cvector = arma::Col<T>;
template<class T>
using rvector = arma::Row<T>;
#endif //NETREG_TYPES_HPP


