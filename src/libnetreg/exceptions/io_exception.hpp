#ifndef NETREG_IOEXCEPTION_HPP
#define NETREG_IOEXCEPTION_HPP

namespace netreg
{
    /**
     * Exception when something when wrong at file OP
     */
    class io_exception : public std::exception
    {

    public:

        io_exception
            (const char *error = "IO-error!")
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
#endif //NETREG_IOEXCEPTION_HPP
