/**
 * Author: Simon Dirmeier
 * Date: 7/31/16
 * Email: simon.dirmeier@bsse.ethz.ch
 */

#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>

int main()
{
    using namespace boost::lambda;
    typedef std::istream_iterator<int> in;

    std::for_each(
        in(std::cin), in(), std::cout << (_1 * 3) << " " );
}