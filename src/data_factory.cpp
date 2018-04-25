/**
 * Author: Simon Dirmeier
 * Date: 25.04.18
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include "data_factory.hpp"

const std::string netreg::data_factory::GAUSSIAN = "gaussian";
const std::string netreg::data_factory::BINOMIAL = "binomial";
const std::string netreg::data_factory::FAMILY_ERROR =
  "Wrong family given. Choose one of " +
  netreg::data_factory::GAUSSIAN + "/" +
  netreg::data_factory::BINOMIAL;
