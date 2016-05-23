#include <CppUTest/TestHarness.h>
#include <vector>
#include "../IndexedVector.hpp"
#include "../Matrix.hpp"
#include "../stats.hpp"
#include "../utility.hpp"

namespace
{
    const double TOL = .001;

    TEST_GROUP(StatsTester)
    {
    };

    TEST(StatsTester, TestVectorSSE)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 2.5;
        const double xval = 1.0;
        std::vector<double> vec(P, yval);
        IndexedVector<double> idx(std::move(vec), -1);
        std::vector<int> testSet = netutil::iota(N, 0);
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        DOUBLES_EQUAL(0.0, netutil::sse(idx, X, Y, testSet, 0), TOL);
    };

    TEST(StatsTester, TestVectorSSE2)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 2.5;
        const double xval = 0.0;
        std::vector<double> vec(P, yval);
        IndexedVector<double> idx(std::move(vec), -1);
        std::vector<int> testSet = netutil::iota(N, 0);
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        const double expected_value_ = yval * yval * N;
        DOUBLES_EQUAL(expected_value_, netutil::sse(idx, X, Y, testSet, 0),
                      TOL);
    };

    TEST(StatsTester, TestMatrixSSE)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 25.0;
        const double xval = 1.0;
        std::vector<int> testSet = netutil::iota(P, 0);
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        Matrix coef(P, Q, yval / Q);
        DOUBLES_EQUAL(0.0, netutil::sse(coef, X, Y, testSet), TOL);
    };

    TEST(StatsTester, TestMatrixSSE2)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 25.0;
        const double xval = 0.0;
        std::vector<int> testSet = netutil::iota(N, 0);
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        Matrix coef(P, Q, yval / Q);
        const double expected_value_ = yval * yval * N * Q;
        DOUBLES_EQUAL(expected_value_, netutil::sse(coef, X, Y, testSet),
                      TOL);
    };

    TEST(StatsTester, TestPartialResidual)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 25.0;
        const double xval = 1.0;
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        std::vector<double> b(P, yval);
        const double expected_value_ = 25 - xval * yval * (P - 1);
        DOUBLES_EQUAL(expected_value_,
                      netutil::partialResidual(X, Y, b, 0, 0, 0),
                      TOL);
    };

    TEST(StatsTester, TestPartialResidual2)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 25.0;
        const double xval = 0.0;
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        std::vector<double> b(P, yval);
        const double expected_value_ = yval;
        DOUBLES_EQUAL(expected_value_,
                      netutil::partialResidual(X, Y, b, 0, 0, 0),
                      TOL);
    };

    TEST(StatsTester, TestPartialResidual3)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 25.0;
        const double xval = 0.0;
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        std::vector<double> b(P);
        for (int i = 0; i < P; ++i)
            b[i] = static_cast<double>(i);
        const int partial_res_idx = 5;
        double expected_value_ = yval;
        for (int j = 0; j < P; ++j)
        {
            if (j != partial_res_idx)
                expected_value_ -= xval * b[j];
        }
        DOUBLES_EQUAL(expected_value_,
                      netutil::partialResidual(X, Y, b, 0, partial_res_idx,
                                               0),
                      TOL);
    };

    TEST(StatsTester, TestMatrixPartialResidual)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 25.0;
        const double xval = 0.0;
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        Matrix B(P, Q, yval);
        const double expected_value_ = 25.0;
        DOUBLES_EQUAL(expected_value_,
                      netutil::partialResidual(X, Y, B, 0, 0, 0),
                      TOL);
    };

    TEST(StatsTester, TestMatrixPartialResidual2)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 25.0;
        const double xval = 1.0;
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        Matrix B(P, Q, yval);
        double expected_value_ = 25.0 - (P - 1) * xval * yval;
        DOUBLES_EQUAL(expected_value_,
                      netutil::partialResidual(X, Y, B, 0, 0, 0),
                      TOL);
    };

    TEST(StatsTester, TestMatrixPartialResidual3)
    {
        const int N = 100;
        const int P = 10;
        const int Q = P;
        const double yval = 25.0;
        const double xval = 1.0;
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        Matrix B(P, Q);
        for (int i = 0; i < P; ++i)
            for (int j = 0; j < Q; ++j)
                B(i, j) = static_cast<double>(i * Q + j);
        const int partial_res_idx = 5;
        const int qi = 0;
        double expected_value_ = 25.0;
        for (int i = 0; i < P; ++i)
            if (i != partial_res_idx)
                expected_value_ -= xval * B(i, 0);
        DOUBLES_EQUAL(expected_value_,
                      netutil::partialResidual(X, Y, B, 0, partial_res_idx,
                                               qi),
                      TOL);
    };

    TEST(StatsTester, TestTransposeResidual4)
    {
        const int N = 100;
        const int P = 10;
        const int Q = 5;
        const double yval = 25.0;
        const double xval = 1.0;
        Matrix Y(N, Q, yval);
        Matrix X(N, P, xval);
        Matrix B(P, Q, yval / P);
        Matrix terr(Q, N);
        netutil::transpose_residual(terr, X, Y, B);
        const double expected_value_ = 0.0;
        for (int i = 0; i < Q; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                DOUBLES_EQUAL(expected_value_, terr(i, j), TOL);
            }
        }
    }

};