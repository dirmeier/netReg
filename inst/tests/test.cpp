/**
 * Author: Simon Dirmeier
 * Date: 7/31/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#define BOOST_TEST_MODULE example
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( free_test_function )
/* Compare with void free_test_function() */
{
        BOOST_TEST( true /* test assertion */ );
}