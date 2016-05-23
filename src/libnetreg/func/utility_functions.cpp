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

    matrix<double> to_matrix(std::vector< std::vector<double> > &elems)
    {
        const int nrow = elems.size();
        const int ncol = elems[0].size();
        matrix<double> m(nrow, ncol);
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