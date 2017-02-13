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

#include "utility_functions.hpp"

#include <cmath>
#include <sstream>

namespace netreg
{
    std::vector<double> to_double(const std::string& str,
                                  const char delim)
    {
        std::vector<double> elems;
        std::stringstream ss(str);
        std::string item;
        while (std::getline(ss, item, delim))
        {
            elems.push_back(atof(item.c_str()));
        }
        return elems;
    }

    std::vector<double> to_double(char *str, const char *spl)
    {
        char *pch;
        pch = strtok(str, spl);
        std::vector<double> v;
        while (pch != NULL)
        {
            v.push_back(atof(pch));
            pch = strtok(NULL, spl);
        }
        return v;
    }

    arma::Mat<double> to_matrix(std::vector< std::vector<double> > &elems)
    {
        const int nrow = elems.size();
        const int ncol = elems[0].size();
        arma::Mat<double> m(nrow, ncol);
        std::vector< std::vector<double> >::iterator row;
        std::vector<double>::iterator col;
        for (row = elems.begin(); row != elems.end(); row++)
        {
            int ri = std::distance(elems.begin(), row);
            for (col = row->begin(); col != row->end(); col++)
            {
                int ci = std::distance(row->begin(), col);
                m(ri, ci) = *col;
            }
        }
        return m;
    }
}
