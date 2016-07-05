#include <CppUTest/TestHarness.h>
#include <vector>
#include "../IndexedVector.hpp"
#include "../Matrix.hpp"
#include "../stats.hpp"
#include "../utility.hpp"
#include "../ccd_util.hpp"

namespace
{
    const double TOL = .001;

    TEST_GROUP(CDDUtilTest)
    {
    };

    TEST(CDDUtilTest, TestUpNorm)
    {
        const int N = 100;
        const int P = 10;
        const double xval = 2.0;
        Matrix X(N, P, xval);
        const double txx = 0.07;
        const int idx = 0;
        std::vector<int> testIdx = netutil::iota(P, 0);
         double val = netutil::up_norm(X, txx, testIdx, idx);
        const double expected_val = txx - P * xval * xval;
        DOUBLES_EQUAL(expected_val, val, TOL);
    };

    TEST(CDDUtilTest, TestPreTXX)
    {
        const int N = 100;
        const int P = 10;
        const double xval = 2.0;
        const int idx = 5;
        const int row = 0;
        Matrix X(N, P, xval);
        Matrix tXX(P, P, 0.0);
        netutil::pre_txx(tXX, X, idx, row);
        const int expected_val = xval * xval;
        for (int i = 0; i < P; ++i)
        {
            for (int j = 0; j < P; ++j)
            {
                if (i == idx)
                {
                        if (j > i) {DOUBLES_EQUAL(0.0, tXX(i, j), TOL);}
                        else                    {DOUBLES_EQUAL(expected_val, tXX(i, j), TOL);}
                }
                else {
                    DOUBLES_EQUAL(0.0, tXX(i, j), TOL);
                }
            }
        }
    };

    TEST(CDDUtilTest, TestPLS)
    {
        const int P = 10;
        const int idx = 3;
        const double txy = 3.0;
        const double val = 2.0;
        Matrix txx(P, P, val);
        const double bval = 2.5;
        std::vector<double> b(P, bval);
        const double expected_value = txy - (P - 1) * val * bval;
        const double pls = netutil::pls(P, idx, txy, txx, b);
        DOUBLES_EQUAL(expected_value, pls, TOL);
    };

    TEST(CDDUtilTest, TestLPLS)
    {
        const int P = 10;
        const int idx = 3;
        const double txy = 3.0;
        const double val = 2.0;
        Matrix txx(P, P, 0.0);
        for (int i = 0; i < P; ++i)
        {
            for (int j = 0; j <= i; ++j)
            {
                txx(i,j) = val;
            }
        }
        const double bval = 2.5;
        std::vector<double> b(P, bval);
        const double expected_value = txy - (P - 1) * val * bval;
        const double pls = netutil::l_pls(P, idx, txy, txx, b);
        DOUBLES_EQUAL(expected_value, pls, TOL);
    };
};
