/**
 * netReg: graph-regularized linear regression models.
 * <p>
 * Copyright (C) 2015 - 2016 Simon Dirmeier
 * <p>
 * This file is part of netReg.
 * <p>
 * netReg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * netReg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with netReg. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author: Simon Dirmeier
 * @email: simon.dirmeier@gmx.de
 */

#include "graph_model_cv_data.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace netreg
{
    void graph_model_cv_data::set_fold_ids()
    {
        if (static_cast<int>(FOLD_IDS_.size()) != DATA_.sample_count())
        {
            FOLD_IDS_.resize(
              static_cast<std::vector<int>::size_type>(DATA_.sample_count()));
        }

        #pragma omp parallel for
        for (int i = 0; i < CVSET_.fold_count(); ++i)
        {
            cv_fold& fold = CVSET_.get_fold(i);
            for (arma::uvec::iterator j = fold.test_set().begin();
                 j != fold.test_set().end(); ++j)
            {
                FOLD_IDS_[*j] = i;
            }
        }
    }
}
