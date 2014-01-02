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

#ifndef SUZERAIN_NDX_HPP
#define SUZERAIN_NDX_HPP

/** @file
 * Symbolic indices, names, and descriptions for 3D Navier--Stokes.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Symbolic indices, identifiers, and descriptions for 3D Navier--Stokes.
 */
namespace ndx {

enum type {

#ifndef SUZERAIN_PARSED_BY_DOXYGEN
    // Force both (e - rho) and (rho - e) to provide usable, signed results.
    ensure_values_are_signed_to_permit_arbitrary_relative_offsets = -1,
#endif

    e   = 0, /**< Signed index for total energy         \f$e   = \rho{}E\f$ */
    mx  = 1, /**< Signed index for streamwise momentum  \f$m_x = \rho{}u\f$ */
    my  = 2, /**< Signed index for wall-normal momentum \f$m_y = \rho{}v\f$ */
    mz  = 3, /**< Signed index for spanwise momentum    \f$m_z = \rho{}w\f$ */
    rho = 4  /**< Signed index for density              \f$\rho         \f$ */
};

/**
 * Short identifiers that permit distinguishing equations.
 * E.g. "rho_E".
 */
extern const array<const char *, ndx::rho + 1u> identifier;

/**
 * Descriptions of the physical quantity associated with each index.
 * E.g. "total energy".
 */
extern const array<const char *, ndx::rho + 1u> description;

/**
 * Compute the index of the (<em>zero-indexed</em>) <tt>i</tt>-th species
 * partial density for <tt>i > 0</tt>.  Because the diluter species partial
 * density if not tracked, <tt>i > 0</tt> is required.  Notice neither
 * \ref identifier nor \ref description is provided for these indices.
 *
 * @tparam Index Type of the index to return.
 * @param  i     Zero-indexed species of interest.
 *
 * @return The equation index for the <tt>i</tt>-th special partial density.
 */
template< typename Index >
Index rho_(const Index i)
{
    assert(i > 0);
    return i + static_cast<Index>(rho);
}

} // namespace ndx

} // namespace suzerain

#endif // SUZERAIN_NDX_HPP
