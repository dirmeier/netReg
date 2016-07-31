/**
 * Author: Simon Dirmeier
 * Date: 7/31/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include "CppUTest/CommandLineTestRunner.h"
#include "CppUTest/TestHarness.h"
int main(int ac, char** av)
{
    return CommandLineTestRunner::RunAllTests(ac, av);
}


TEST_GROUP(FirstTestGroup)
    {
    };

TEST(FirstTestGroup, FirstTest)
{
FAIL("Fail me!");
}

TEST(FirstTestGroup, SecondTest)
{
    STRCMP_EQUAL("hello", "world");
    LONGS_EQUAL(1, 2);
    CHECK(false);
}

