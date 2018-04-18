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

#ifndef NETREG_EDGENET_MODEL_SELECTION_HPP
#define NETREG_EDGENET_MODEL_SELECTION_HPP

#include <map>

#include "params.hpp"
#include "graph_model_cv_data.hpp"
#include "cross_validation.hpp"

namespace netreg
{
    /**
     * Find the set of optimal shrinkage parameters for a edge-penalized
     * regression
     * model. Set is calculated using cross-validation.
     *
     * @template validator a validator
     * @template deviance the deviance of some distribution
     *
     * @param data an object that holds all required data for the model
     * @param pars a parameter object containing all relevant informations
     *
     * @returns returns a map of shrinkage parameters
     */
    template<template<typename ...> class validator>
    std::map<std::string, double> model_selection(
      graph_model_cv_data& data, params& pars);
}

#endif
