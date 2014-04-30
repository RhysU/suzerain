//--------------------------------------------------------------------------
//
// Copyright (C) 2012-2014 Rhys Ulerich
// Copyright (C) 2012-2014 The PECOS Development Team
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

#ifndef SUZERAIN_SPECIFICATION_LARGO_HPP
#define SUZERAIN_SPECIFICATION_LARGO_HPP

/** @file
 * Provides \ref specification_largo.
 */

#include <suzerain/common.hpp>
#include <suzerain/baseflow.hpp>
#include <suzerain/largo_formulation.hpp>

struct largo_workspace;

namespace suzerain {

/**
 * Holds parameters defining Largo-based slow growth problems.
 */
class specification_largo
{
public:

    /** Construct an instance with the given values. */
    explicit specification_largo(
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

    /** Growth rate of amplitudes for the mean. */
    std::vector<real_t> gramp_mean;

    /** Growth rate of amplitudes for the rms. */
    std::vector<real_t> gramp_rms;

    /**
     * Should slow growth computations ignore all fluctuations?
     * Meant as a debugging tool emphatically \e not preserved across restarts.
     */
    bool ignore_fluctuations;

    /**
     * Should slow growth computations consider #gramp_mean to be zero?
     * Meant as a debugging tool emphatically \e not preserved across restarts.
     */
    bool ignore_gramp_mean;

    /**
     * Should slow growth computations consider #gramp_rms to be zero?
     * Meant as a debugging tool emphatically \e not preserved across restarts.
     */
    bool ignore_gramp_rms;

};

} // namespace suzerain

#endif // SUZERAIN_SPECIFICATION_LARGO_HPP
