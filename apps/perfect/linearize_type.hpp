//--------------------------------------------------------------------------
//
// Copyright (C) 2011-2014 Rhys Ulerich
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

#ifndef SUZERAIN_PERFECT_LINEARIZE_TYPE_HPP
#define SUZERAIN_PERFECT_LINEARIZE_TYPE_HPP

/** @file
 * Declarations related to linearization type selection.
 */

namespace suzerain {

namespace perfect {

/** Provides scoping semantics for linearize::type */
namespace linearize {

/**
 * What type of hybrid implicit/explicit linearization is employed?
 * All valid types will evaluate to true in a boolean context.
 */
enum type {
    rhome_xyz = 1, ///< Density, momentum, and total energy across X, Y, and Z
    rhome_y,       ///< Density, momentum, and total energy implicit in Y only
    none           ///< No linearization implying a fully explicit treatment
};

} // namespace linearize

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_LINEARIZE_TYPE_HPP */
