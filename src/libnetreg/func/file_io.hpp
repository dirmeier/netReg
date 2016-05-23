/**
 * Author: Simon Dirmeier
 * Date: 3/17/16
 * Email: uccd_.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_FILE_IO_HPP
#define NETREG_FILE_IO_HPP

#include <string>
#include <vector>

#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#include <armadillo>

#include "../struc/types.hpp"
#include "../data/modeldata/linear_model_data.hpp"
#include "../struc/optim/pareto_optimal_point.hpp"

namespace netreg
{
    matrix<double> read_matrix(const std::string &file);

    std::vector <std::vector<double>>
        read_content(const std::string &file);

    void write_results(linear_model_data &data, std::string &output);

    void write_results(pareto_optimal_point<std::string, double> &data,
                       std::string &output);
};
#endif //NETREG_FILE_IO_HPP
