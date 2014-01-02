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

#ifndef SUZERAIN_BASEFLOW_HPP
#define SUZERAIN_BASEFLOW_HPP

/** @file
 * Provides \ref baseflow_interface and implementations.
 */

#include <suzerain/common.hpp>
#include <suzerain/largo_state.hpp>

// Forward declarations
struct suzerain_radialflow_solution;

namespace suzerain {

// Forward declarations
class bspline;

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
     * @param[out] base   Flow state at \f$y\f$ stored per \ref largo_state.
     * @param[out] dybase Wall-normal derivatives of flow state at \f$y\f$.
     *                    Stored per \ref largo_state.
     * @param[out] dxbase Streamwise derivatives of flow state at \f$y\f$.
     *                    Stored per \ref largo_state.
     */
    virtual void conserved(const real_t      y,
                           real_t *       base,
                           real_t *     dybase,
                           real_t *     dxbase) const = 0;

    /**
     * Compute baseflow pressure and its spatial derivatives at some position.
     *
     * @param[in ] y   Wall-normal coordinate \f$y\f$.
     * @param[out] P   Base flow pressure at \f$y\f$.
     * @param[out] dyP Wall-normal derivative of pressure at \f$y\f$.
     * @param[out] dxP Streamwise derivative of pressure at \f$y\f$.
     */
    virtual void pressure(const real_t    y,
                          real_t &        P,
                          real_t &      dyP,
                          real_t &      dxP) const = 0;

};

/**
 * Provides slow-growth-ready baseflow information for uniform flows.
 * That is, flows for which all spatial derivatives are zero.
 */
class baseflow_uniform : public virtual baseflow_interface
{
public:
    /**
     * Construct an instance
     *
     * Member #x must be initialized prior to using
     * conserved() or pressure().
     */
    baseflow_uniform();

    void conserved(const real_t      y,
                   real_t *       base,
                   real_t *     dybase,
                   real_t *     dxbase) const;

    void pressure(const real_t    y,
                  real_t &        P,
                  real_t &      dyP,
                  real_t &      dxP) const;

    /**
     * Each variable in the baseflow is represented by one entry in \c x.
     * Pressure information is stored in the final column.
     */
    VectorXr x;
};

/**
 * Provides slow-growth-ready baseflow information from polynomial fits.
 */
class baseflow_polynomial : public virtual baseflow_interface
{
public:
    /**
     * Construct an instance
     *
     * Members #x and #dx must be initialized prior to using
     * conserved() or pressure().
     */
    baseflow_polynomial();

    void conserved(const real_t      y,
                   real_t *       base,
                   real_t *     dybase,
                   real_t *     dxbase) const;

    void pressure(const real_t    y,
                  real_t &        P,
                  real_t &      dyP,
                  real_t &      dxP) const;

    /**
     * Polynomial coefficients for the fields.  Each variable in the baseflow
     * is represented by one column in \c x.  Each column encodes one
     * polynomial in the format the <a
     * href="http://www.gnu.org/software/gsl/manual/html_node/Polynomial-Evaluation.html">
     * GSL Polynomial Evaluation</a> routines expect.  Derivatives in the
     * \f$y\f$ direction are taken from the same coefficient representation.
     * Pressure information is stored in the final column.  This member must
     * have the same number of columns as #dx but the number of rows may
     * differ.
     */
    MatrixXXr x;

    /**
     * Polynomial coefficients for the field \f$x\$ derivatives.  These are
     * stored analogously to #x.  In particular, this member must have the same
     * number of columns as #x but the number of rows may differ.
     */
    MatrixXXr dx;
};

/**
 * Provides map-based lookup of arbitrary baseflow profiles established at \e
 * runtime.  This may be used with, e.g., \ref suzerain_radialflow_solver
 * and \ref suzerain_radialflow_cartesian_conserved to compute baseflow
 * information on the collocation points for a given \ref suzerain::bspline.
 */
class baseflow_map : public virtual baseflow_interface
{
public:

    baseflow_map();

    /**
     * Look up conserved state results from \c table.
     * If not found and possible, interpolate for the results.
     * \throw std::out_of_range on requests requiring extrapolation.
     */
    void conserved(const real_t      y,
                   real_t *       base,
                   real_t *     dybase,
                   real_t *     dxbase) const;

    /**
     * Look up pressure results from \c table.
     * If not found and possible, interpolate for the results.
     * \throw std::out_of_range on requests requiring extrapolation.
     */
    void pressure(const real_t    y,
                  real_t &        P,
                  real_t &      dyP,
                  real_t &      dxP) const;

    /** Each entry in the table consists of this information. */
    struct row { largo_state base, dybase, dxbase; };

    /** Keys are \c y coordinates returning \c row instances. */
    typedef std::map<real_t, row> table_type;

    /** Tracks known data with keys being \c y coordinates. */
    table_type table;
};

} // namespace suzerain

#endif // SUZERAIN_BASEFLOW_HPP
