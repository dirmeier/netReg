#include <CppUTest/TestHarness.h>
#include <vector>
#include <algorithm>
#include "../utility.hpp"

namespace
{
    TEST_GROUP(UtilityTester)
    {
    };

    TEST(UtilityTester, TestIota)
    {
        const int le = 10;
        const int start = 5;
        std::vector<int> vec = netutil::iota(le, start);
        CHECK_EQUAL(le, vec.size());
    };

    TEST(UtilityTester, TestIota2)
    {
        const int le = 10;
        const int start = 5;
        std::vector<int> vec = netutil::iota(le, start);
        int runner = 5;
        for (int i = 0; i < le; ++i)
        {
            CHECK_EQUAL(vec[i], runner);
            runner++;
        }
    };

    TEST(UtilityTester, TestShuffle)
    {
        const int le = 10;
        const int start = 5;
        std::vector<int> vec = netutil::iota(le, start);
        netutil::shuffle(vec);
        int smaller_count = 0;
        for (std::vector<int>::size_type i = 0; i < vec.size() - 1; i++)
        {
            smaller_count += vec[i] <= vec[i + 1] ? 1 : 0;
        }
        CHECK_TRUE(smaller_count < le - 1);
    };

    TEST(UtilityTester, TestShuffle2)
    {
        const int le = 10;
        const int start = 5;
        std::vector<int> vec = netutil::shuffle(le, start);
        int smaller_count = 0;
        for (std::vector<int>::size_type i = 0; i < vec.size() - 1; i++)
        {
            smaller_count += vec[i] <= vec[i + 1] ? 1 : 0;
        }
        CHECK_TRUE(smaller_count < le - 1);
    };

    TEST(UtilityTester, TestShuffle3)
    {
        const int le = 10;
        const int start = 5;
        std::vector<int> vec = netutil::shuffle(le, start);
        std::sort(vec.begin(), vec.end());
        int runner = 5;
        for (int i = 0; i < le; ++i)
        {
            CHECK_EQUAL(vec[i], runner);
            runner++;
        }
    };
};
