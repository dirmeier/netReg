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


#ifndef NETREG_PARAMS_HPP
#define NETREG_PARAMS_HPP

namespace netreg
{
    class params
    {
    public:
        params& lambda(double lambda)
        {
            lambda_ = lambda;
            return *this;
        }

        params& psigx(double psigx)
        {
            psigx_ = psigx;
            return *this;
        }

        params& psigy(double psigy)
        {
            psigy_ = psigy;
            return *this;
        }

        params& do_lambda(bool do_lambda)
        {
            do_lambda_ = do_lambda;
            return *this;
        }

        params& do_psigx(bool do_psigx)
        {
            do_psigx_ = do_psigx;
            return *this;
        }

        params& do_psigy(bool do_psigy)
        {
            do_psigy_ = do_psigy;
            return *this;
        }

        params& thresh(double thresh)
        {
            thresh_ = thresh;
            return *this;
        }

        params& niter(int niter)
        {
            niter_ = niter;
            return *this;
        }

        params& optim_niter(int optim_niter)
        {
            optim_niter_ = optim_niter;
            return *this;
        }

        params& optim_epsilon(double optim_epsilon)
        {
            optim_epsilon_ = optim_epsilon;
            return *this;
        }

    private:
        double lambda_;
        double psigx_;
        double psigy_;

        bool do_lambda_;
        bool do_psigx_;
        bool do_psigy_;

        double thresh_;
        int niter_;

        int optim_niter_;
        double optim_epsilon_;
    };
}

#endif  // NETREG_PARAMS
