//--------------------------------------------------------------------------
//
// Copyright (C) 2010-2014 Rhys Ulerich
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

#ifndef SUZERAIN_PERFECT_OPERATOR_COMMON_BLOCK_HPP
#define SUZERAIN_PERFECT_OPERATOR_COMMON_BLOCK_HPP

/** @file
 * Implements \ref operator_common_block.
 */

#include <suzerain/common.hpp>

#include "implicits.hpp"
#include "instantaneous.hpp"
#include "linearize_type.hpp"
#include "references.hpp"

namespace suzerain {

namespace perfect {

/**
 * Storage for holding quantities computed during nonlinear operator
 * application and implicit constraint application which either are required
 * for linear operator application or for statistics sampling purposes.
 */
class operator_common_block
    : public virtual boost::noncopyable
{
public:

    /** Default constructor.  Use \ref set_zero to resize prior to use. */
    operator_common_block();

    /**
     * Determines the extent of the implicit treatment
     * by the paired linear and nonlinear operators.
     */
    linearize::type linearization;

    /** Call \c set_zero on each of \c ref, \c sub, and \c imp. */
    void set_zero(int Ny);

    /**
     * Reference quantities to be updated by the nonlinear operator
     * at the beginning of each Runge--Kutta substep.
     */
    references ref;

    /**
     * Instantaneous quantities to be updated by the nonlinear operator
     * at the beginning of each Runge--Kutta substep.
     */
    instantaneous sub;

    /**
     * The implicit quantities, stored as collocation point values in \ref
     * imp and updated substep-by-substep by \ref constraint::treatment.
     */
    implicits imp;

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_OPERATOR_COMMON_BLOCK_HPP */
