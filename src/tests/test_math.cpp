#include <CppUTest/TestHarness.h>
#include <vector>
#include "../math.hpp"

namespace
{
    const double TOL_ = .001;

    TEST_GROUP(MathTester)
    {
    };

    TEST(MathTester, TestConverges)
    {
        const int le = 10;
        std::vector<double> vec(le, 1.0);
        std::vector<double> vec2(le, 1.0);
        bool converges_ = netutil::converges(vec, vec2, TOL_);
        CHECK_TRUE(converges_);
    };

    TEST(MathTester, TestConverges2)
    {
        const int le = 10;
        std::vector<double> vec(le, 1.0);
        std::vector<double> vec2(le, 1.1);
        bool converges_ = netutil::converges(vec, vec2, TOL_);
        CHECK_FALSE(converges_);
    };

    TEST(MathTester, TestConverges6)
    {
        const int le = 10;
        const double xval = 0.0;
        const double yval = TOL_ / le;
        std::vector<double> vec(le, xval);
        std::vector<double> vec2(le, yval);
        bool converges_ = netutil::converges(vec, vec2, TOL_);
        CHECK_FALSE(converges_);
    };

    TEST(MathTester, TestConverges7)
    {
        const int le = 10;
        const double xval = 0.0;
        const double yval = TOL_ / le - 0.000001;
        std::vector<double> vec(le, xval);
        std::vector<double> vec2(le, yval);
        bool converges_ = netutil::converges(vec, vec2, TOL_);
        CHECK_TRUE(converges_);
    };

    TEST(MathTester, TestConverges3)
    {
        const int N = 10;
        const int M = 10;
        const double val = 0.1;
        Matrix X(N, M, val);
        Matrix Y(N, M, val);
        bool converges_ = netutil::converges(X, Y, TOL_);
        CHECK_TRUE(converges_);
    };

    TEST(MathTester, TestConverges4)
    {
        const int N = 10;
        const int M = 10;
        const double xval = 0.0;
        const double yval = TOL_ / (N * M);
        Matrix X(N, M, xval);
        Matrix Y(N, M, yval);
        bool converges_ = netutil::converges(X, Y, TOL_);
        CHECK_FALSE(converges_);
    };

    TEST(MathTester, TestConverges5)
    {
        const int N = 10;
        const int M = 10;
        const double xval = 0.0;
        const double yval = TOL_ / (N * M) - 0.000001;
        Matrix X(N, M, xval);
        Matrix Y(N, M, yval);
        bool converges_ = netutil::converges(X, Y, TOL_);
        CHECK_TRUE(converges_);
    };

    TEST(MathTester, TestConverges8)
    {
        const int N = 10;
        const int M = 10;
        const double xval = 0.0;
        const double yval = (TOL_ / M) - 0.00001;
        Matrix X(N, M, xval);
        Matrix Y(N, M, yval);
        for (int i = 0; i < M; ++i)
        {
            bool converges_ = netutil::converges(X, Y, i, TOL_);
            CHECK_TRUE(converges_);
        }
    };

    TEST(MathTester, TestConverges9)
    {
        const int N = 10;
        const int M = 10;
        const double xval = 0.0;
        const double yval = (TOL_ / M);
        Matrix X(N, M, xval);
        Matrix Y(N, M, yval);
        for (int i = 0; i < M; ++i)
        {
            bool converges_ = netutil::converges(X, Y, i, TOL_);
            CHECK_FALSE(converges_);
        }
    };

    TEST(MathTester, TestSoft)
    {
        const double s = 0;
        const double lalph = 1.0;
        const double norm = 1.0;
        const double expected_val_ = 0.0;
        DOUBLES_EQUAL(expected_val_, netutil::softnorm(s, lalph, norm),
                      TOL_);
    }

    TEST(MathTester, TestSoft1)
    {
        const double s = 0.11;
        const double lalph = 0.1;
        const double norm = 1.0;
        const double expected_val_ = s - lalph;
        DOUBLES_EQUAL(expected_val_, netutil::softnorm(s, lalph, norm),
                      TOL_);
    }

    TEST(MathTester, TestSoft2)
    {
        const double s = -0.11;
        const double lalph = 0.1;
        const double norm = 1.0;
        const double expected_val_ = s + lalph;
        DOUBLES_EQUAL(expected_val_, netutil::softnorm(s, lalph, norm),
                      TOL_);
    }

    TEST(MathTester, TestSoft3)
    {
        const double s = 0.1;
        const double lalph = 0.1;
        const double norm = 2.0;
        const double expected_val_ = 0.0;
        DOUBLES_EQUAL(expected_val_, netutil::softnorm(s, lalph, norm),
                      TOL_);
    }

    TEST(MathTester, TestSoft4)
    {
        const double s = 0.1;
        const double lalph = 0.01;
        const double norm = 2.0;
        const double expected_val_ = (s - lalph) / norm;
        DOUBLES_EQUAL(expected_val_, netutil::softnorm(s, lalph, norm),
                      TOL_);
    }

    TEST(MathTester, TestAbsSum)
    {
        const int N = 10;
        const int M = 10;
        const double xval = 0.0;
        const double yval = 1.0;
        Matrix X(N, M, xval);
        Matrix Y(N, M, yval);
        const double sum = netutil::abs_sum(X, Y);
        const double expected_val_ = yval * M * N;
        DOUBLES_EQUAL(expected_val_, sum, TOL_);
    };

    TEST(MathTester, TestAbsSum2)
    {
        const int N = 10;
        const int M = 10;
        const double xval = 1.0;
        const double yval = 1.0;
        Matrix X(N, M, xval);
        Matrix Y(N, M, yval);
        const double sum = netutil::abs_sum(X, Y);
        const double expected_val_ = 0.0;
        DOUBLES_EQUAL(expected_val_, sum, TOL_);
    };

    TEST(MathTester, TestAbsSum3)
    {
        const int N = 10;
        const int M = 10;
        const double xval = 1.0;
        const double yval = 1.0;
        Matrix X(N, M, xval);
        Matrix Y(N, M, yval);
        for (int i = 0; i < M; ++i)
        {
            const double sum = netutil::abs_sum(X, Y, i);
            const double expected_val_ = 0.0;
            DOUBLES_EQUAL(expected_val_, sum, TOL_);
        }
    };

    TEST(MathTester, TestAbsSum4)
    {
        const int N = 10;
        const int M = 10;
        const double xval = 0.0;
        const double yval = 1.0;
        Matrix X(N, M, xval);
        Matrix Y(N, M, yval);
        for (int i = 0; i < M; ++i)
        {
            const double sum = netutil::abs_sum(X, Y, i);
            const double expected_val_ = yval * N;
            DOUBLES_EQUAL(expected_val_, sum, TOL_);
        }
    };

    TEST(MathTester, TestAbsSum5)
    {
        const int N = 10;
        const double xval = 0.0;
        const double yval = 1.0;
        std::vector<double> v1(N, xval);
        std::vector<double> v2(N, yval);
        const double sum = netutil::abs_sum(v1, v2);
        const double expected_val_ = yval * N;
        DOUBLES_EQUAL(expected_val_, sum, TOL_);

    };

};
