//--------------------------------------------------------------------------
//
// Copyright (C) 2012 Rhys Ulerich
// Copyright (C) 2012 The PECOS Development Team
// Please see http://pecos.ices.utexas.edu for more information on PECOS.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
 *
 * The last index \e must be <tt>rho</tt>!  In the multi-species cases, mixture
 * density should be in \ref rho and the <tt>i</tt>-th species partial density
 * in <tt>rho + (i-1)</tt>.  That is, \ref rho is the "diluter" species.  Many
 * routines will assume the equations are stored in this order when attempting
 * to walk memory linearly.
 */
namespace ndx {

// Anonymous enum to declare our state variable storage indices.
enum {
    e,    /**< Index for storing total energy         \f$e   = \rho{}E\f$ */
    mx,   /**< Index for storing streamwise momentum  \f$m_x = \rho{}u\f$ */
    my,   /**< Index for storing wall-normal momentum \f$m_y = \rho{}v\f$ */
    mz,   /**< Index for storing spanwise momentum    \f$m_z = \rho{}w\f$ */
    rho   /**< Index for storing density              \f$      \rho   \f$ */
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

} // namespace ndx

} // namespace suzerain

#endif // SUZERAIN_NDX_HPP
