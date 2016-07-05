/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "graph_functions.hpp"

#include <vector>
#include <cmath>

namespace netreg
{
    matrix<double> laplacian(matrix<double> &m)
    {
        const int N_ROW_ = static_cast<int>(m.n_rows);
        const int N_COLS_ = static_cast<int>(m.n_cols);
        matrix<double> lapla(N_ROW_, N_COLS_);
        // node degrees of source matrix
        std::vector<double> degrees(N_ROW_);
        for (int i = 0; i < N_ROW_; ++i)
        {
            double rowSum = 0.0;
            for (int j = 0; j < N_COLS_; ++j)
                rowSum += m(i, j);
            degrees[i] = rowSum;
        }
        // calculate normalized Laplacian matrix of source
        for (int i = 0; i < N_ROW_; ++i)
        {
            for (int j = 0; j < N_COLS_; ++j)
            {
                if (i == j && degrees[i] != 0)
                    lapla(i, j) = 1 - (m(i, j) / degrees[i]);
                else if (i != j && m(i, j) != 0)
                    lapla(i, j) =
                        -m(i, j) / std::sqrt(degrees[i] * degrees[j]);
                else
                    lapla(i, j) = 0.0;
            }
        }
        return lapla;
    }
}