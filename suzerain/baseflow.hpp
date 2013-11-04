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

#ifndef SUZERAIN_BASEFLOW_HPP
#define SUZERAIN_BASEFLOW_HPP

/** @file
 * Provides \ref baseflow_interface and implementations.
 */

#include <suzerain/common.hpp>

namespace suzerain {

/**
 * Abstract interface for providing slow-growth-ready baseflow information.
 */
class baseflow_interface
{
public:

    baseflow_interface();

    /** Virtual destructor for interface-like class */
    virtual ~baseflow_interface();

    /**
     * Compute baseflow state and its spatial derivatives at some position.
     *
     * @param[in ] y      Wall-normal coordinate \f$y\f$.
     * @param[out] base   Base flow state at \f$y\f$ stored per \ref largo_state.
     * @param[out] dybase Wall-normal derivatives of flow state at \f$y\f$.
     *                    Stored per \ref largo_state.
     * @param[out] dxbase Streamwise derivatives of flow state at \f$y\f$.
     *                    Stored per \ref largo_state.
     */
    virtual void get_baseflow(const real_t    y,
                              real_t *     base,
                              real_t *   dybase,
                              real_t *   dxbase) const = 0;

    /**
     * Compute baseflow pressure and its spatial derivatives at some position.
     *
     * @param[in ] y      Wall-normal coordinate \f$y\f$.
     * @param[out] base   Base flow pressure at \f$y\f$.
     * @param[out] dybase Wall-normal derivative of pressure at \f$y\f$.
     * @param[out] dxbase Streamwise derivative of pressure at \f$y\f$.
     */
    virtual void get_baseflow_pressure(const real_t    y,
                                       real_t &    Pbase,
                                       real_t &  dyPbase,
                                       real_t &  dxPbase) const = 0;

};

} // namespace suzerain

#endif // SUZERAIN_BASEFLOW_HPP
