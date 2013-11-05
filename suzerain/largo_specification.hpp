//--------------------------------------------------------------------------
//
// Copyright (C) 2012, 2013 Rhys Ulerich
// Copyright (C) 2012, 2013 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This file is part of Suzerain.
//
// Suzerain is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Suzerain is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Suzerain.  If not, see <http://www.gnu.org/licenses/>.
//
//--------------------------------------------------------------------------

#ifndef SUZERAIN_LARGO_SPECIFICATION_HPP
#define SUZERAIN_LARGO_SPECIFICATION_HPP

/** @file
 * Provides \ref largo_specification.
 */

#include <suzerain/common.hpp>
#include <suzerain/baseflow.hpp>
#include <suzerain/largo_formulation.hpp>

struct largo_workspace;

namespace suzerain {

/**
 * Holds parameters defining Largo-based slow growth problems.
 */
class largo_specification
{
public:

    /** Construct an instance with the given values. */
    explicit largo_specification(
            largo_formulation formulation = largo_formulation::disable,
            real_t grdelta                = std::numeric_limits<real_t>::quiet_NaN(),
            largo_workspace * workspace   = NULL);

    /** Which \ref largo_formulation is in use? */
    largo_formulation formulation;

    /** Growth rate of reference thickness \f$\Delta\f$. */
    real_t grdelta;

    /** Pointer to largo workspace. */
    largo_workspace * workspace;

    /** Baseflow information. */
    shared_ptr<baseflow_interface> baseflow;

};

} // namespace suzerain

#endif // SUZERAIN_LARGO_SPECIFICATION_HPP
