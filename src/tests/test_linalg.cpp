#include <CppUTest/TestHarness.h>
#include <vector>
#include "../Matrix.hpp"
#include "../linalg.hpp"
#include "../stats.hpp"

namespace
{
    const double TOL_ = .001;

    TEST_GROUP(LinAlgTester)
    {
    };

    TEST(LinAlgTester, TestTranspose)
    {
        const int N = 10;
        const int M = 5;
        Matrix X(N, M);
        Matrix Y(M, N);
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                X(i, j) = i * M + j;
            }
        }
        netutil::transpose(Y, X);
        for (int i = 0; i < M; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                DOUBLES_EQUAL(Y(i, j), X(j, i), TOL_);
            }
        }
    };

    TEST(LinAlgTester, TestMatMul)
    {
        const int N = 10;
        const int M = 5;
        Matrix X(N, M, 1.0);
        Matrix Y(M, N, 1.0);
        netutil::transpose(Y, X);
        Matrix Z(N, N);
        netutil::matmul(Z, X, Y);
        const double expected_value = 5.0;
        for (int i = 0; i < M; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                DOUBLES_EQUAL(expected_value, Z(i, j), TOL_);
            }
        }
    };

    TEST(LinAlgTester, TestMatMul2)
    {
        const int N = 10;
        const int M = 0;
        Matrix X(N, M, 1.0);
        Matrix Y(M, N, 1.0);
        netutil::transpose(Y, X);
        Matrix Z(N, N);
        netutil::matmul(Z, X, Y);
        const double expected_value = 0.0;
        for (int i = 0; i < M; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                DOUBLES_EQUAL(expected_value, Z(i, j), TOL_);
            }
        }
    };

    TEST(LinAlgTester, TestMatMul3)
    {
        const int N = 10;
        const int M = 5;
        const double xval = .2;
        const double yval = .1;
        Matrix X(N, M, xval);
        Matrix Y(M, N, yval);
        Matrix Z(N, N);
        netutil::matmul(Z, X, Y);
        const double expected_value = xval * yval * M;
        for (int i = 0; i < M; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                DOUBLES_EQUAL(expected_value, Z(i, j), TOL_);
            }
        }
    };

    TEST(LinAlgTester, TestAbsMaxDotProduct)
    {
        const int N = 10;
        const int M = 5;
        const double smval = .1;
        const double bgval = .2;
        Matrix X(N, M, smval);
        Matrix Y(N, M, smval);
        for (int k = 0; k < N; ++k)
            X(k, 0) = bgval;
        for (int k = 0; k < M; ++k)
            Y(0, k) = bgval;
        const double expected_value = (0.2 * 0.2) + (0.2 * 0.1) * (N - 1);
        const double dot_product = netutil::maxAbsDotProduct(X, Y);
        DOUBLES_EQUAL(expected_value, dot_product, TOL_);
    };

    TEST(LinAlgTester, TestAbsDotProduct)
    {
        const int N = 10;
        const int M = 5;
        const double smval = .1;
        const double bgval = .2;
        Matrix X(N, M, smval);
        Matrix Y(N, M, smval);
        for (int k = 0; k < N; ++k)
            X(k, 0) = bgval;
        for (int k = 0; k < M; ++k)
            Y(0, k) = bgval;
        for (int i = 0; i < M; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                double expected_value = 0.0;
                for (int k = 0; k < N; ++k)
                    expected_value += X(k, i) * Y(k, j);
                const double dot_product = netutil::absDotproduct(X, Y, i, j);
                DOUBLES_EQUAL(expected_value, dot_product, TOL_);
            }
        }
    };

    TEST(LinAlgTester, TestAbsDotProduct1)
    {
        const int N = 10;
        const int M = 5;
        const double smval = .1;
        const double bgval = .2;
        Matrix X(N, M, smval);
        Matrix Y(N, M, smval);
        for (int k = 0; k < N; ++k)
            X(k, 0) = bgval;
        for (int k = 0; k < M; ++k)
            Y(0, k) = bgval;
        for (int i = 0; i < M; ++i)
        {
            for (int j = 0; j < M; ++j)
            {
                double expected_value = 0.0;
                for (int k = 0; k < N; ++k)
                    expected_value += X(k, i) * Y(k, j);
                const double dot_product = netutil::absDotproduct(X, Y, i, j);
                DOUBLES_EQUAL(expected_value, dot_product, TOL_);
            }
        }
    };

    TEST(LinAlgTester, TestAbsColMaxDotProduct)
    {
        const int N = 10;
        const int M = 5;
        const double smval = .1;
        const double bgval = .2;
        Matrix X(N, M, smval);
        Matrix Y(N, M, smval);
        for (int k = 0; k < N; ++k)
            X(k, 0) = bgval;
        double expected_value = .2 * .1 * N;
        for (int j = 0; j < M; ++j)
        {
            const double dot_product = netutil::maxColWiseAbsDotProduct(X, Y, j);
            DOUBLES_EQUAL(expected_value, dot_product, TOL_);
        }

    };

    TEST(LinAlgTester, TestIntercept)
    {
        const int N = 10;
        const int P = 5;
        const int Q = 3;
        const double smval = 0;
        const double bgval = 1;
        Matrix X(N, P, smval);
        Matrix Y(N, Q, bgval);
        Matrix B(P, Q, bgval / P);
        Matrix terr(Q, N);
        netutil::transpose_residual(terr, X, Y, B);
        std::vector<double> intr(Q);
        netutil::intercept(intr, terr);
        const double expected_value_ = 1;
        for (int i = 0; i < Q; ++i)
        {
            DOUBLES_EQUAL(expected_value_, intr[i], TOL_);
        }
    };

    TEST(LinAlgTester, TestIntercept1)
    {
        const int N = 10;
        const int P = 5;
        const int Q = 3;
        const double smval = 0;
        const double bgval = 1;
        Matrix X(N, P, smval);
        Matrix Y(N, Q, bgval);
        Matrix B(P, Q, bgval / P);
        std::vector<double> intr(Q);
        netutil::intercept(intr, X, Y, B);
        const double expected_value_ = 1;
        for (int i = 0; i < Q; ++i)
        {
            DOUBLES_EQUAL(expected_value_, intr[i], TOL_);
        }
    };

};
