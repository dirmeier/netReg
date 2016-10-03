/**
 * netReg: graph-regularized linear regression models.
 * <p>
 * Copyright (C) 2015 - 2016 Simon Dirmeier
 * <p>
 * This file is part of netReg.
 * <p>
 * netReg is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * netReg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with netReg. If not, see <http://www.gnu.org/licenses/>.
 *
 * @author: Simon Dirmeier
 * @email: simon.dirmeier@gmx.de
 */
#ifndef NETREG_ILLEGAL_ARGUMENT_EXCEPTION_HPP
#define NETREG_ILLEGAL_ARGUMENT_EXCEPTION_HPP

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
