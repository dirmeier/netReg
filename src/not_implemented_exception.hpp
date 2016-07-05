/**
 * Author: Simon Dirmeier
 * Email: netreg@simon-dirmeier.net
 */

#ifndef NETREG_NOTIMPLEMENTEDEXCEPTION_HPP
#define NETREG_NOTIMPLEMENTEDEXCEPTION_HPP

#include <exception>
#include <string>

namespace netreg
{
    /**
     * Exception when a function is called that is not yet implemented.
     */
    class not_implemented_exception : public std::exception
    {

    public:
        not_implemented_exception
            (const char *error = "Functionality not yet implemented!")
        {
            errorMessage = error;
        }

        const char *what() const noexcept
        {
            return errorMessage.c_str();
        }

    private:

        std::string errorMessage;
    };
}
#endif //NETREG_NOTIMPLEMENTEDEXCEPTION_HPP
