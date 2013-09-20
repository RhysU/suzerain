//--------------------------------------------------------------------------
//
// Copyright (C) 2013 Rhys Ulerich
// Copyright (C) 2013 The PECOS Development Team
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

#ifndef SUZERAIN_PERFECT_LARGO_STATE_HPP
#define SUZERAIN_PERFECT_LARGO_STATE_HPP

/** @file
 * Implementation of nonlinear Navier--Stokes spatial operators.
 */

#include <suzerain/common.hpp>

namespace suzerain {

namespace perfect {

/**
 * A Largo-related convenience class for manipulating double-valued buffers
 * containing \f$\rho\f$, \f$\rho{}u\f$, \f$\rho{}v\f$, \f$\rho{}w\f$,
 * \f$\rho{}E\f$, and \f$p\f$ in that order.  Largo-based slow growth
 * computations require conserved state packed in this fashion.  Pressure
 * \f$p\f$ is wholly redundant but Largo requires it for certain slow growth
 * forcing models.
 */
union largo_state
{
public:

    /** Initialize with all zeros such that \ref trivial() holds. */
    largo_state() { std::memset(this, 0, sizeof(largo_state)); }

    /** Initialize with scalars per order of suzerain:ndx::type. */
    largo_state(double e, double mx, double my, double mz, double rho, double p)
        : rho(rho), mx(mx), my(my), mz(mz), e(e), p(p) {}

    // Storage directly accessible through public members
    struct {
        double rho;   /**< Density              \f$\rho   \f$ */
        double mx;    /**< Streamwise momentum  \f$\rho{}u\f$ */
        double my;    /**< Wall-normal momentum \f$\rho{}v\f$ */
        double mz;    /**< Spanwise momentum    \f$\rho{}w\f$ */
        double e;     /**< Total energy         \f$\rho{}E\f$ */
        double p;     /**< Pressure             \f$p      \f$ */
    };

private:

    /** Access to the packed data gated through rescale() or as_is() below. */
    double state[6];

public:

    // Helpers for computing commonly-accessed results
    double u() const { return mx / rho; } /**< Computes \f$u = m_x/\rho \f$ */
    double v() const { return my / rho; } /**< Computes \f$v = m_y/\rho \f$ */
    double w() const { return mz / rho; } /**< Computes \f$w = m_z/\rho \f$ */
    double E() const { return e  / rho; } /**< Computes \f$E = e  /\rho \f$ */

    /** Is the state uniformly zero? */
    bool trivial() const
    {
        using namespace std;
#pragma warning(push,disable:1572)
        return 0 == abs(rho) + abs(mx) + abs(my) + abs(mz) + abs(e) + abs(p);
#pragma warning(pop)
    }

    /**
     * Largo is assumed to work correctly when the Euler equations are
     * nondimensionalized using \f$l_0\f$, \f$u_0 = a_0\f$, and \f$p_0 =
     * \rho_0 u_u^2\f$ as that produces nondimensional equations with the
     * same form as the dimensional ones.
     *
     * When \f$ a_0 \ne u_0 \f$ one has non-unit \f$\mbox{Ma}\f$ and both
     * total energy and pressure must be scaled by \f$\mbox{Ma}^{-2}\f$
     * when calling Largo and scaled by \f$\mbox{Ma}^2\f$ on return.
     *
     * The following is a helper class for performing that conversion.  The
     * constructor scales <tt>s->e</tt> and <tt>s->p</tt> by some constant \c C
     * and the destructor undoes that scaling factor.  The implicit conversion
     * to <tt>double*</tt> permits passing \c helper instances directly to
     * Largo API calls.  This crazy concoction being useful relies heavily on
     * C++ temporary object lifetime semantics.
     */
    class rescaler
    {
    public:
        rescaler (largo_state *s, double C)
            : s(s), C(C)
        { s->e *= C; s->p *= C; }

        ~rescaler()
        { s->e /= C; s->p /= C; }

        operator double*()
        { return s->state; }

    private:
        largo_state *s;
        double C;
    };

    friend class rescaler;

    /**
     * Use this method when state needs to be passed to/from Largo.
     * @see Discussion at #rescaler.
     */
    rescaler rescale(double inv_Ma2) { return rescaler(this, inv_Ma2); }

    /** Use this method when state needs does not go to/from Largo. */
    double * as_is() { return this->state; }

};

} // namespace perfect

} // namespace suzerain

#endif  /* SUZERAIN_PERFECT_LARGO_STATE_HPP */
