/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */
#ifndef NETREG_PARETOOPTIMALPOINT_HPP
#define NETREG_PARETOOPTIMALPOINT_HPP

#include <map>
#include <vector>
#include <utility>

namespace netreg
{
    /**
     * Representation of a pareto optimal point (for the sake of using this term)
     *
     * @template T class of intern name representations (e.g. std::string)
     * @template U class of the point data-types (e.g. double)
     */
    template<typename T, typename U>
    class pareto_optimal_point
    {
    public:
        /**
         * Add a pareto optimal value to the set.
         */
        void put(T key, U value)
        {
            keys_.push_back(key);
            values_.push_back(value);
        }

        /**
         * Getter for the name-value pair of a parameter.
         *
         * @return name-value pair
         */
        std::pair <T, U> &operator[](const int idx)
        {
            entry_.first = keys_[idx];
            entry_.second = values_[idx];
            return entry_;
        }

        /**
         * Getter for the number of parameters of the point.
         *
         * @return the number of parameters
         */
        int npar()
        {
            return values_.size();
        }

    private:
        std::vector <T> keys_;   // names of the values
        std::vector <U> values_; // minimal values
        std::pair <T, U> entry_; // a single name-value pair instantiated in order to avoid constant reallocation.
    };

}

#endif //NETREG_PARETOOPTIMALPOINT_HPP
