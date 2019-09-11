/**
 * Author: Simon Dirmeier
 * Date: 02/10/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#ifndef NETREG_FAMILY_HPP
#define NETREG_FAMILY_HPP

#include <cstdint>

namespace netreg
{
    enum class family : std::int8_t
    {
        BINOMIAL = 0,
        GAUSSIAN = 1,
        NONE     = -1
    };
}

#endif  // NETREG_FAMILY_HPP
