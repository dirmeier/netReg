#ifndef NETREG_ILLEGALARGUMENTEXCEPTION_HPP
#define NETREG_ILLEGALARGUMENTEXCEPTION_HPP

#include <exception>
#include <string>

namespace netreg
{
    /**
     * Exception when wrong function parameters where provided
     */
    class illegal_argument_exception : public std::exception
    {
    public:

        illegal_argument_exception
            (const char *error = "Wrong arguments given")
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
#endif //NETREG_ILLEGALARGUMENTEXCEPTION_HPP
