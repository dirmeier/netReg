/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#include "graph_functions.hpp"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include <vector>
#include <cmath>

namespace netreg
{
    matrix<double> laplacian(const double * x, const int n, const int m)
    {
        std::vector<double> degrees(n);
        #pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            double rowSum = 0.0;
            for (int j = 0; j < m; ++j)
                rowSum += x[i + n * j];
            degrees[i] = rowSum;
        }
        double *  laplacian = new double[n * m];
        // calculate normalized Laplacian matrix of source
        #pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                if (i == j && degrees[i] != 0)
                    laplacian[i + n * j] = 1 - (x[i + n * j] / degrees[i]);
                else if (i != j && x[i + n * j] != 0)
                    laplacian[i + n * j] =
                        -x[i + n * j] / std::sqrt(degrees[i] * degrees[j]);
                else
                    laplacian[i + n * j] = 0.0;
            }
        }
        matrix<double> lap(laplacian, n, m, true, true);
        delete [] laplacian;
        return lap;
    }
}